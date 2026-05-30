## -----------------------------------------------------------------------------
#| include: false
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = has_pkg
)


## -----------------------------------------------------------------------------
#| label: setup
library(TemporalHazard)
library(survival)
library(ggplot2)


## -----------------------------------------------------------------------------
#| label: explore-data
data(avc)
avc <- na.omit(avc)

# Candidate covariates
candidates <- c("age", "status", "mal", "com_iv", "inc_surg", "orifice")


## -----------------------------------------------------------------------------
#| label: explore-logit
logit_fits <- lapply(candidates, function(var) {
  fmla <- as.formula(paste("dead ~", var))
  glm(fmla, data = avc, family = binomial())
})
names(logit_fits) <- candidates

# Which covariates are significant univariably?
p_values <- vapply(logit_fits, function(f) {
  coef(summary(f))[2, "Pr(>|z|)"]
}, numeric(1))

data.frame(covariate = names(p_values),
           p_value   = round(p_values, 4),
           row.names = NULL)


## -----------------------------------------------------------------------------
#| label: fig-explore-age
#| fig-cap: "Exploratory: age vs. mortality with LOESS smooth"
#| fig-width: 7
#| fig-height: 4
ggplot(avc, aes(age, dead)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, colour = "#0072B2",
              linewidth = 0.8) +
  labs(x = "Age at operation (months)", y = "Death (0/1)") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: calibrate-age
age_cal <- hzr_calibrate(
  x      = avc$age,
  event  = avc$dead,
  groups = 10,
  link   = "logit"
)
print(age_cal)


## -----------------------------------------------------------------------------
#| label: km-baseline
km <- hzr_kaplan(time = avc$int_dead, status = avc$dead, conf_level = 0.95)
print(km)


## -----------------------------------------------------------------------------
#| label: nelson-baseline
nel <- hzr_nelson(time = avc$int_dead, event = avc$dead, conf_level = 0.95)
print(nel)


## -----------------------------------------------------------------------------
#| label: fit-base
# Fit the base model that bootstrapping will use
fit <- hazard(
  Surv(int_dead, dead) ~ age + status + mal + com_iv,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.20, nu = 1.0, rep(0, 4)),
  fit   = TRUE,
  control = list(maxit = 500)
)

# Prediction grid at a median-risk profile
t_grid <- seq(0.01, max(avc$int_dead) * 0.9, length.out = 100)
base_nd <- data.frame(time = t_grid, age = median(avc$age),
                      status = 2, mal = 0, com_iv = 0)

surv_point <- predict(fit, newdata = base_nd, type = "survival")


## -----------------------------------------------------------------------------
#| label: bootstrap
set.seed(42)
boot <- hzr_bootstrap(fit, n_boot = 30)  # kept small for vignette build time
print(boot)


## -----------------------------------------------------------------------------
#| label: gof-summary
gof <- hzr_gof(fit)
print(gof)


## -----------------------------------------------------------------------------
#| label: deciles
deciles <- hzr_deciles(fit, time = 120, groups = 10)
print(deciles)


## -----------------------------------------------------------------------------
#| label: fig-deciles
#| fig-cap: "Observed vs. expected mortality by decile of risk"
#| fig-width: 7
#| fig-height: 4.5
decile_df <- as.data.frame(deciles)
cal_long <- rbind(
  data.frame(group = decile_df$group,
             rate  = decile_df$observed_rate,
             Series = "Observed"),
  data.frame(group = decile_df$group,
             rate  = decile_df$expected_rate,
             Series = "Expected")
)
ggplot(cal_long, aes(group, rate, fill = Series)) +
  geom_col(position = position_dodge(width = 0.7), alpha = 0.8) +
  scale_fill_manual(values = c("Observed" = "#56B4E9",
                                "Expected" = "#E69F00")) +
  scale_x_continuous(breaks = seq_len(nrow(decile_df))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Risk decile (1 = lowest)", y = "Mortality rate", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: sensitivity
# Define reference and high-risk profiles
ref_profile <- data.frame(
  age    = median(avc$age),
  status = 2,
  mal    = 0,
  com_iv = 0
)

high_risk <- data.frame(
  age    = quantile(avc$age, 0.90),
  status = 4,
  mal    = 1,
  com_iv = 1
)


## -----------------------------------------------------------------------------
#| label: fig-sensitivity
#| fig-cap: "Sensitivity analysis: reference vs. high-risk profile"
#| fig-width: 7
#| fig-height: 4.5
sens_curves <- do.call(rbind, lapply(
  list("Reference" = ref_profile, "High risk" = high_risk),
  function(p) {
    nd <- p[rep(1, length(t_grid)), ]
    nd$time <- t_grid
    data.frame(time = t_grid,
               survival = predict(fit, newdata = nd, type = "survival") * 100,
               Profile = deparse(substitute(p)))
  }
))
# Fix profile labels
sens_curves$Profile <- rep(c("Reference", "High risk"),
                           each = length(t_grid))
sens_curves$Profile <- factor(sens_curves$Profile,
                              levels = c("Reference", "High risk"))

ggplot(sens_curves, aes(time, survival, colour = Profile)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("Reference" = "#009E73",
                                 "High risk" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Months after AVC repair", y = "Freedom from death (%)",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

