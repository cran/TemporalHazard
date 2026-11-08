## -----------------------------------------------------------------------------
#| include: false
# R CMD build installs the package into R_BUILD_TEMPLIB. A build launched
# through callr (rcmdcheck, devtools) passes Quarto's R session a profile
# that resets .libPaths() without it, and every chunk skipped (#276).
templib <- Sys.getenv("R_BUILD_TEMPLIB")
if (dir.exists(templib)) .libPaths(c(templib, .libPaths()))
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)
# On CI a package that cannot load must fail the build, not skip every
# chunk and pass over nothing (#276). loadNamespace() says why it failed.
if (nzchar(Sys.getenv("CI"))) for (p in c("TemporalHazard", "ggplot2")) loadNamespace(p)
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
#| label: kul-data
data(cabgkul)
str(cabgkul)


## -----------------------------------------------------------------------------
#| label: kul-weibull
fit_kul <- hazard(
  Surv(int_dead, dead) ~ 1,
  data  = cabgkul,
  dist  = "weibull",
  theta = c(mu = 0.10, nu = 1.0),
  fit   = TRUE
)

fit_kul


## -----------------------------------------------------------------------------
#| label: fig-kul-km
#| fig-cap: "Weibull parametric survival vs. Kaplan-Meier (CABG, KU Leuven)"
#| fig-width: 7
#| fig-height: 4.5
t_grid <- seq(0.01, max(cabgkul$int_dead) * 0.9, length.out = 200)
nd     <- data.frame(time = t_grid)
surv   <- predict(fit_kul, newdata = nd, type = "survival") * 100

km     <- survfit(Surv(int_dead, dead) ~ 1, data = cabgkul)
km_df  <- data.frame(time = km$time, survival = km$surv * 100)

ggplot() +
  geom_step(data = km_df, aes(time, survival, colour = "Kaplan-Meier"),
            linewidth = 0.5) +
  geom_line(data = data.frame(time = t_grid, survival = surv),
            aes(time, survival, colour = "Weibull"), linewidth = 1) +
  scale_colour_manual(values = c("Weibull" = "#0072B2",
                                 "Kaplan-Meier" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Months after CABG", y = "Freedom from death (%)",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: avc-data
data(avc)
avc <- na.omit(avc)
str(avc)


## -----------------------------------------------------------------------------
#| label: avc-mv
fit_avc <- hazard(
  Surv(int_dead, dead) ~ age + status + mal + com_iv + inc_surg + orifice,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.20, nu = 1.0, rep(0, 6)),
  fit   = TRUE,
  control = list(maxit = 500)
)

fit_avc


## -----------------------------------------------------------------------------
#| label: avc-mp
fit_mp <- hazard(
  Surv(int_dead, dead) ~ 1,
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
#| label: fig-avc-compare
#| fig-cap: "Single-phase Weibull vs. multiphase model against Kaplan-Meier (AVC)"
#| fig-width: 7
#| fig-height: 4.5
t_grid <- seq(0.01, max(avc$int_dead) * 0.95, length.out = 200)
nd     <- data.frame(time = t_grid)

km_avc <- survfit(Surv(int_dead, dead) ~ 1, data = avc)
km_df  <- data.frame(time = km_avc$time, survival = km_avc$surv * 100)

fit_wb <- hazard(
  Surv(int_dead, dead) ~ 1, data = avc, dist = "weibull",
  theta = c(mu = 0.20, nu = 1.0), fit = TRUE
)

surv_wb <- predict(fit_wb, newdata = nd, type = "survival") * 100
surv_mp <- predict(fit_mp, newdata = nd, type = "survival") * 100

plot_df <- rbind(
  data.frame(time = t_grid, survival = surv_wb, Model = "Single Weibull"),
  data.frame(time = t_grid, survival = surv_mp, Model = "Multiphase (2-phase)")
)

ggplot() +
  geom_step(data = km_df, aes(time, survival), colour = "grey50",
            linewidth = 0.5) +
  geom_line(data = plot_df, aes(time, survival, colour = Model),
            linewidth = 1) +
  scale_colour_manual(values = c("Single Weibull" = "#E69F00",
                                 "Multiphase (2-phase)" = "#0072B2")) +
  scale_y_continuous(limits = c(0, 100)) +
  annotate("text", x = max(t_grid) * 0.6, y = 95, label = "KM (grey)",
           size = 3, colour = "grey50") +
  labs(x = "Months after AVC repair", y = "Freedom from death (%)",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: valves-death
data(valves)
valves <- na.omit(valves)

fit_death <- hazard(
  Surv(int_dead, dead) ~ age_cop + nyha + mechvalv,
  data  = valves,
  dist  = "weibull",
  theta = c(mu = 0.10, nu = 1.0, rep(0, 3)),
  fit   = TRUE,
  control = list(maxit = 500)
)

fit_death


## -----------------------------------------------------------------------------
#| label: valves-pve
fit_pve <- hazard(
  Surv(int_pve, pve) ~ age_cop + nve + mechvalv,
  data  = valves,
  dist  = "weibull",
  theta = c(mu = 0.02, nu = 1.0, rep(0, 3)),
  fit   = TRUE,
  control = list(maxit = 500)
)

fit_pve


## -----------------------------------------------------------------------------
#| label: interval-simulate
set.seed(101)
n         <- 400
visit_gap <- 6   # months between visits

# True event times: Weibull mu = 0.025, nu = 1.0 (median ~28 months)
true_time <- rexp(n, rate = 0.025)  # nu=1 → exponential

n_visits  <- sample(4:8, n, replace = TRUE)

status     <- integer(n)
time       <- numeric(n)
time_lower <- numeric(n)
time_upper <- numeric(n)

for (i in seq_len(n)) {
  visits <- seq(visit_gap, n_visits[i] * visit_gap, by = visit_gap)
  t_ev   <- true_time[i]
  brk    <- findInterval(t_ev, c(0, visits))  # which interval contains t_ev

  if (brk > length(visits)) {
    # Event after last visit: right-censored.
    # time_lower = 0: with H(time_lower) = H(0) = 0, the right-censored
    # contribution reduces to -(H(time) - H(0)) = -H(time), which is correct.
    status[i]     <- 0L
    time[i]       <- max(visits)
    time_lower[i] <- 0
    time_upper[i] <- max(visits)
  } else {
    # Event within the visit window: interval-censored in (lower, upper)
    lower_bound   <- if (brk == 1L) 0 else visits[brk - 1L]
    upper_bound   <- visits[brk]
    status[i]     <- 2L
    time[i]       <- lower_bound
    time_lower[i] <- lower_bound
    time_upper[i] <- upper_bound
  }
}

table(status)


## -----------------------------------------------------------------------------
#| label: interval-fit
fit_ic <- hazard(
  time       = time,
  status     = status,
  time_lower = time_lower,
  time_upper = time_upper,
  dist  = "weibull",
  theta = c(mu = 0.025, nu = 1.0),
  fit   = TRUE
)

fit_ic


## -----------------------------------------------------------------------------
#| label: interval-vs-naive
# Naive: treat each interval-censored row as an exact event at time_upper
# (the visit when the death was first recorded)
fit_naive <- hazard(
  time   = ifelse(status == 2L, time_upper, time),
  status = ifelse(status == 2L, 1L,         status),
  dist  = "weibull",
  theta = c(mu = 0.025, nu = 1.0),
  fit   = TRUE
)

# Compare MLEs; true parameters are mu = 0.025, nu = 1.0
rbind(
  interval_censored = round(coef(fit_ic)[1:2],   4),
  naive_exact_upper = round(coef(fit_naive)[1:2], 4),
  truth             = c(mu = 0.025, nu = 1.0)
)


## -----------------------------------------------------------------------------
#| label: starting-values
nel <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)

# Drop the zero-time boundary point before log-transforming
nel_clean <- nel[nel$time > 0 & nel$cumhaz > 0, ]
log_t <- log(nel_clean$time)
log_H <- log(nel_clean$cumhaz)

lm_fit  <- lm(log_H ~ log_t)
nu_hat  <- unname(coef(lm_fit)[2])
mu_hat  <- exp(unname(coef(lm_fit)[1]) / nu_hat)
c(mu = round(mu_hat, 4), nu = round(nu_hat, 4))


## -----------------------------------------------------------------------------
#| label: fig-cumhaz-shape
#| fig-cap: "Nelson-Aalen cumulative hazard for CABGKUL: three-phase shape visible as two kinks"
#| fig-width: 7
#| fig-height: 4
ggplot(data.frame(time = nel$time, cumhaz = nel$cumhaz),
       aes(time, cumhaz)) +
  geom_step(colour = "#D55E00", linewidth = 0.7) +
  labs(x = "Months after CABG", y = "Cumulative hazard") +
  theme_minimal()


## -----------------------------------------------------------------------------
#| label: overparameterized
fit_3ph <- hazard(
  Surv(int_dead, dead) ~ 1,
  data   = avc,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant"),
    late     = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                          fixed = "shapes")
  ),
  fit = TRUE, control = list(n_starts = 5, maxit = 1000)
)

# Scale magnitudes. A tiny exp(log_mu) does not by itself mark an unneeded
# phase; judge each phase by its contribution, from predict(..., decompose = TRUE)
log_mu_idx <- grep("log_mu", names(coef(fit_3ph)))
round(exp(coef(fit_3ph)[log_mu_idx]), 6)


## -----------------------------------------------------------------------------
#| label: nstarts
fit_robust <- hazard(
  Surv(int_dead, dead) ~ 1,
  data   = avc,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant")
  ),
  fit     = TRUE,
  control = list(n_starts = 10, maxit = 1000)
)
fit_robust$fit$converged

