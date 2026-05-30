## -----------------------------------------------------------------------------
#| include: false
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  eval = has_pkg
)


## -----------------------------------------------------------------------------
#| label: setup
library(TemporalHazard)
if (has_ggplot2) library(ggplot2)

data(avc)
avc <- na.omit(avc)
str(avc)


## -----------------------------------------------------------------------------
#| label: km-life-table
km <- hzr_kaplan(time = avc$int_dead, status = avc$dead)
head(km)


## -----------------------------------------------------------------------------
#| label: fig-km-baseline
#| fig-cap: "Kaplan-Meier survival estimate for AVC patients (n = 310)"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
km_df <- data.frame(
  time     = km$time,
  survival = km$survival * 100,
  lower    = km$cl_lower * 100,
  upper    = km$cl_upper * 100
)

ggplot(km_df, aes(time, survival)) +
  geom_step(colour = "#D55E00", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              stat = "identity", alpha = 0.15, fill = "#D55E00") +
  labs(x = "Months after AVC repair", y = "Survival (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: weibull-fit
fit_weib <- hazard(
  survival::Surv(int_dead, dead) ~ 1,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.05, nu = 0.5),
  fit   = TRUE
)
summary(fit_weib)


## -----------------------------------------------------------------------------
#| label: fig-weibull-overlay
#| fig-cap: "Single Weibull vs. Kaplan-Meier"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
t_grid <- seq(0.01, max(avc$int_dead) * 0.9, length.out = 200)
surv_weib <- predict(fit_weib,
                     newdata = data.frame(time = t_grid),
                     type = "survival") * 100

ggplot() +
  geom_step(data = km_df, aes(time, survival),
            colour = "#D55E00", linewidth = 0.6) +
  geom_line(data = data.frame(time = t_grid, survival = surv_weib),
            aes(time, survival), colour = "#0072B2", linewidth = 0.8) +
  labs(x = "Months after AVC repair", y = "Survival (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: multiphase-fit
fit_mp <- hazard(
  survival::Surv(int_dead, dead) ~ 1,
  data   = avc,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant")
  ),
  fit     = TRUE,
  control = list(n_starts = 5, maxit = 1000)
)
summary(fit_mp)


## -----------------------------------------------------------------------------
#| label: fig-multiphase-overlay
#| fig-cap: "Two-phase parametric model vs. Kaplan-Meier"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
surv_mp <- predict(fit_mp,
                   newdata = data.frame(time = t_grid),
                   type = "survival") * 100

ggplot() +
  geom_step(data = km_df, aes(time, survival),
            colour = "#D55E00", linewidth = 0.6) +
  geom_line(data = data.frame(time = t_grid, survival = surv_mp),
            aes(time, survival), colour = "#0072B2", linewidth = 0.8) +
  labs(x = "Months after AVC repair", y = "Survival (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: fig-decomposed
#| fig-cap: "Phase decomposition of cumulative hazard"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
decomp <- predict(fit_mp,
                  newdata = data.frame(time = t_grid),
                  type = "cumulative_hazard",
                  decompose = TRUE)

decomp_df <- data.frame(
  time  = t_grid,
  total = decomp[, "total"],
  early = decomp[, "early"],
  const = decomp[, "constant"]
)

ggplot(decomp_df, aes(x = time)) +
  geom_line(aes(y = total), linewidth = 0.9, colour = "black") +
  geom_line(aes(y = early), linewidth = 0.7, colour = "#E69F00",
            linetype = "dashed") +
  geom_line(aes(y = const), linewidth = 0.7, colour = "#56B4E9",
            linetype = "dashed") +
  labs(x = "Months after AVC repair", y = "Cumulative hazard") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: screening
covariates <- c("age", "status", "mal", "com_iv", "orifice",
                "inc_surg", "opmos")

screen_results <- data.frame(
  variable = covariates,
  coef     = numeric(length(covariates)),
  p_value  = numeric(length(covariates))
)

for (i in seq_along(covariates)) {
  v <- covariates[i]
  fml <- as.formula(paste("dead ~", v))
  lg <- glm(fml, data = avc, family = binomial)
  s <- summary(lg)$coefficients
  if (nrow(s) >= 2) {
    screen_results$coef[i] <- s[2, 1]
    screen_results$p_value[i] <- s[2, 4]
  }
}

screen_results <- screen_results[order(screen_results$p_value), ]
screen_results$significant <- ifelse(screen_results$p_value < 0.05,
                                      "*", "")
screen_results


## -----------------------------------------------------------------------------
#| label: age-calibrate
cal_age <- hzr_calibrate(x = avc$age, event = avc$dead,
                          groups = 10, link = "logit")
cal_age


## -----------------------------------------------------------------------------
#| label: fig-calibrate-age
#| fig-cap: "Decile calibration: age vs. logit(P(death))"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
ggplot(cal_age, aes(mean, link_value)) +
  geom_point(size = 3, colour = "#0072B2") +
  geom_line(colour = "#0072B2", alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#D55E00", linetype = "dashed") +
  labs(x = "Age at repair (months, decile midpoint)",
       y = "logit(P(death))") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: multivariable-fit
fit_mv <- hazard(
  survival::Surv(int_dead, dead) ~ age + status + mal + com_iv,
  data   = avc,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant")
  ),
  fit     = TRUE,
  control = list(n_starts = 5, maxit = 1000)
)
summary(fit_mv)


## -----------------------------------------------------------------------------
#| label: stepwise-forward
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE)
base_mp <- hazard(
  survival::Surv(int_dead, dead) ~ 1,
  data   = avc,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant")
  ),
  fit     = TRUE,
  control = list(n_starts = 3, maxit = 500)
)

fit_step <- hzr_stepwise(
  base_mp,
  scope     = list(
    early    = ~ age + status + mal + com_iv,
    constant = ~ age + status + mal + com_iv
  ),
  data      = avc,
  direction = "both",
  criterion = "wald",
  trace     = FALSE,
  control   = list(n_starts = 2, maxit = 500)
)
fit_step


## -----------------------------------------------------------------------------
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE)
fit_step$steps[, c("step_num", "action", "variable", "phase",
                   "p_value", "aic")]


## -----------------------------------------------------------------------------
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE)
logLik_manual <- fit_mv$fit$objective
logLik_step   <- fit_step$fit$objective
c(manual = logLik_manual, stepwise = logLik_step,
  aic_manual = 2 * length(fit_mv$fit$theta) - 2 * logLik_manual,
  aic_step   = 2 * length(fit_step$fit$theta) - 2 * logLik_step)


## -----------------------------------------------------------------------------
#| label: fig-prediction
#| fig-cap: "Reference patient survival with 95% delta-method CI and KM overlay"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
ref_patient <- data.frame(
  time   = t_grid,
  age    = median(avc$age),
  status = 1,
  mal    = 0,
  com_iv = 0
)

surv_ref <- predict(fit_mv, newdata = ref_patient,
                    type = "survival", se.fit = TRUE)

ref_df <- data.frame(
  time     = t_grid,
  survival = surv_ref$fit * 100,
  lower    = surv_ref$lower * 100,
  upper    = surv_ref$upper * 100
)

ggplot() +
  geom_step(data = km_df, aes(time, survival),
            colour = "#D55E00", linewidth = 0.6, alpha = 0.7) +
  geom_ribbon(data = ref_df, aes(time, ymin = lower, ymax = upper),
              fill = "#0072B2", alpha = 0.15) +
  geom_line(data = ref_df, aes(time, survival),
            colour = "#0072B2", linewidth = 0.8) +
  labs(x = "Months after AVC repair", y = "Survival (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: fig-sensitivity
#| fig-cap: "Survival by risk profile: low-risk (blue) vs. high-risk (red)"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
low_risk <- data.frame(
  time   = t_grid,
  age    = quantile(avc$age, 0.25),
  status = 1,
  mal    = 0,
  com_iv = 0
)

high_risk <- data.frame(
  time   = t_grid,
  age    = quantile(avc$age, 0.90),
  status = 4,
  mal    = 1,
  com_iv = 1
)

surv_lo <- predict(fit_mv, newdata = low_risk,
                   type = "survival", se.fit = TRUE)
surv_hi <- predict(fit_mv, newdata = high_risk,
                   type = "survival", se.fit = TRUE)

sens_df <- data.frame(
  time     = rep(t_grid, 2),
  survival = c(surv_lo$fit, surv_hi$fit) * 100,
  lower    = c(surv_lo$lower, surv_hi$lower) * 100,
  upper    = c(surv_lo$upper, surv_hi$upper) * 100,
  profile  = rep(c("Low risk", "High risk"), each = length(t_grid))
)

ggplot(sens_df, aes(time, survival, colour = profile, fill = profile)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Low risk" = "#0072B2",
                                  "High risk" = "#D55E00")) +
  scale_fill_manual(values = c("Low risk" = "#0072B2",
                                "High risk" = "#D55E00")) +
  labs(x = "Months after AVC repair", y = "Survival (%)",
       colour = NULL, fill = NULL) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: fig-deciles
#| fig-cap: "Observed vs. expected event rates by risk decile"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
cal <- hzr_deciles(fit_mv, time = max(avc$int_dead))
print(cal)


## -----------------------------------------------------------------------------
#| label: fig-decile-plot
#| fig-cap: "Decile calibration: observed (bars) vs. expected (points) event rates"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
ggplot(cal, aes(x = group)) +
  geom_col(aes(y = observed_rate), fill = "#56B4E9", alpha = 0.7) +
  geom_point(aes(y = expected_rate), colour = "#D55E00", size = 3) +
  geom_line(aes(y = expected_rate), colour = "#D55E00", linewidth = 0.5) +
  scale_x_continuous(breaks = seq_len(nrow(cal))) +
  labs(x = "Risk decile (1 = lowest)", y = "Event rate",
       caption = paste("Overall chi-sq =",
                        round(attr(cal, "overall")$chi_sq, 2),
                        ", p =",
                        format.pval(attr(cal, "overall")$p_value,
                                     digits = 3))) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: gof-summary
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE)
gof <- hzr_gof(fit_mv)
print(gof)


## -----------------------------------------------------------------------------
#| label: fig-gof-survival
#| fig-cap: "Parametric (blue) vs. Kaplan-Meier (orange) survival"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
ggplot(gof, aes(x = time)) +
  geom_step(aes(y = km_surv * 100), colour = "#D55E00", linewidth = 0.6) +
  geom_line(aes(y = par_surv * 100), colour = "#0072B2", linewidth = 0.8) +
  labs(x = "Months after AVC repair", y = "Survival (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: fig-gof-events
#| fig-cap: "Cumulative observed (orange) vs. expected (blue) events"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
ggplot(gof, aes(x = time)) +
  geom_line(aes(y = cum_observed), colour = "#D55E00", linewidth = 0.8) +
  geom_line(aes(y = cum_expected), colour = "#0072B2",
            linewidth = 0.8, linetype = "dashed") +
  labs(x = "Months after AVC repair", y = "Cumulative events") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: fig-gof-residual
#| fig-cap: "Conservation-of-events residual (expected minus observed) over time"
#| eval: !expr requireNamespace("TemporalHazard", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
ggplot(gof, aes(x = time, y = residual)) +
  geom_line(linewidth = 0.7, colour = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = "Months after AVC repair",
       y = "Residual (expected - observed)") +
  theme_minimal()

