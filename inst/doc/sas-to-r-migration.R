## -----------------------------------------------------------------------------
#| include: false
# R CMD build installs the package into R_BUILD_TEMPLIB. A build launched
# through callr (rcmdcheck, devtools) passes Quarto's R session a profile
# that resets .libPaths() without it, and every chunk skipped (#276).
templib <- Sys.getenv("R_BUILD_TEMPLIB")
if (dir.exists(templib)) .libPaths(c(templib, .libPaths()))
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE)
# On CI a package that cannot load must fail the build, not skip every
# chunk and pass over nothing (#276). loadNamespace() says why it failed.
if (nzchar(Sys.getenv("CI"))) loadNamespace("TemporalHazard")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = has_pkg
)


## -----------------------------------------------------------------------------
#| label: setup
library(TemporalHazard)


## -----------------------------------------------------------------------------
#| label: mapping-table
knitr::kable(
  TemporalHazard::hzr_argument_mapping(),
  caption = "Formal argument map: SAS HAZARD/C → hazard()",
  col.names = c(
    "SAS Statement", "Legacy Input", "C Concept",
    "R Parameter", "Required", "Expected Type",
    "Transform Rule", "Status", "Notes"
  )
)


## -----------------------------------------------------------------------------
#| eval: false
# fit <- hazard(..., data = avcs)


## -----------------------------------------------------------------------------
#| eval: false
# fit <- hazard(
#   time   = avcs$INT_DEAD,
#   ...
# )


## -----------------------------------------------------------------------------
#| eval: false
# fit <- hazard(
#   status = avcs$DEAD,   # 1 = event, 0 = censored
#   ...
# )


## -----------------------------------------------------------------------------
#| eval: false
# fit <- hazard(
#   ...,
#   dist   = "multiphase",
#   phases = list(
#     early    = hzr_phase("cdf", t_half = 0.1905077, nu = 1.437416, m = 1,
#                          fixed = "m"),       # FIXM: hold M at its start
#     constant = hzr_phase("constant")
#   ),
#   theta = c(log(0.3504743),                  # MUE
#             log(0.1905077), 1.437416, 1,     # THALF, NU, M
#             log(4.391673e-07))               # MUC
# )


## -----------------------------------------------------------------------------
#| eval: false
# fit <- hazard(
#   Surv(int_dead, dead) ~ 1,
#   data   = avcs,
#   dist   = "multiphase",
#   phases = list(
#     early    = hzr_phase("cdf", t_half = 0.1905077, nu = 1.437416, m = 1,
#                          formula = ~ age + com_iv + mal + opmos + op_age +
#                            status,
#                          fixed = "m"),
#     constant = hzr_phase("constant", formula = ~ inc_surg + orifice + status)
#   ),
#   theta = c(
#     log(0.3504743), log(0.1905077), 1.437416, 1,  # MUE, THALF, NU, M
#     -0.03205774, 1.336675, 0.6872028,             # EARLY: AGE, COM_IV, MAL,
#     -0.01963377, 0.0002086689, 0.5169533,         #   OPMOS, OP_AGE, STATUS
#     log(4.391673e-07),                            # MUC
#     1.375285, 3.11765, 1.054988                   # CONSTANT: INC_SURG,
#   ),                                              #   ORIFICE, STATUS
#   fit = TRUE
# )


## -----------------------------------------------------------------------------
#| label: macro-demos-setup
library(survival)
data(avc)
avc <- na.omit(avc)

# Base parametric fit for downstream diagnostics
fit <- hazard(
  Surv(int_dead, dead) ~ age + status + mal + com_iv,
  data = avc, dist = "weibull",
  theta = c(mu = 0.2, nu = 1, rep(0, 4)),
  fit = TRUE, control = list(maxit = 500)
)


## -----------------------------------------------------------------------------
#| label: kaplan-demo
km <- hzr_kaplan(time = avc$int_dead, status = avc$dead)
head(km, 4)


## -----------------------------------------------------------------------------
#| label: nelson-demo
nel <- hzr_nelson(time = avc$int_dead, event = avc$dead)
head(nel, 4)


## -----------------------------------------------------------------------------
#| label: gof-demo
head(hzr_gof(fit), 4)


## -----------------------------------------------------------------------------
#| label: deciles-demo
hzr_deciles(fit, time = 120, groups = 10)


## -----------------------------------------------------------------------------
#| label: calibrate-demo
hzr_calibrate(avc$age, avc$dead, groups = 10, link = "logit")


## -----------------------------------------------------------------------------
#| label: bootstrap-demo
set.seed(1)
hzr_bootstrap(fit, n_boot = 20)  # small for vignette build


## -----------------------------------------------------------------------------
#| label: competing-demo
data(valves)
# Synthesize a competing-risks indicator for illustration.
valves$ev <- with(valves, ifelse(dead == 1, 1L,
                                   ifelse(pve == 1,  2L,
                                   ifelse(reop == 1, 3L, 0L))))
head(hzr_competing_risks(valves$int_dead, valves$ev), 4)


## -----------------------------------------------------------------------------
#| label: avc-example
#| warning: false
data(avc)   # reload: the diagnostics above dropped the incomplete rows
avcs <- avc

# PROC STANDARD REPLACE: fill missing INC_SURG with the column mean
avcs$inc_surg[is.na(avcs$inc_surg)] <- mean(avcs$inc_surg, na.rm = TRUE)

fit <- hazard(
  Surv(int_dead, dead) ~ 1,                        # TIME, EVENT
  data   = avcs,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.1905077, nu = 1.437416, m = 1,
                         formula = ~ age + com_iv + mal + opmos + op_age +
                           status,                 # EARLY
                         fixed = "m"),             # FIXM
    constant = hzr_phase("constant",
                         formula = ~ inc_surg + orifice + status)  # CONSTANT
  ),
  theta = c(                                       # PARMS, EARLY, CONSTANT
    log(0.3504743), log(0.1905077), 1.437416, 1,
    -0.03205774, 1.336675, 0.6872028, -0.01963377, 0.0002086689, 0.5169533,
    log(4.391673e-07),
    1.375285, 3.11765, 1.054988
  ),
  control = list(conserve = TRUE),                 # CONSERVE
  fit = TRUE
)

summary(fit)


## -----------------------------------------------------------------------------
#| eval: false
#| label: predict-example
# # Linear predictor (X %*% theta)
# eta <- predict(fit, type = "linear_predictor")
# 
# # Hazard scale (exp(eta))
# hz  <- predict(fit, type = "hazard")


## -----------------------------------------------------------------------------
#| eval: false
#| label: translate-example
# job <- hzr_translate_sas(
#   system.file("extdata", "hz-example.sas", package = "TemporalHazard"),
#   out_dir = "translated"
# )


## -----------------------------------------------------------------------------
#| eval: false
#| label: translate-output
# # chunk: data
# if (!exists("AVCS")) {
#   stop("This job read ", "AVCS", " from a SAS DATA step, which ",
#        "hzr_translate_sas() does not translate. Assign ", "AVCS",
#        " before rendering.")
# }
# 
# # chunk: status
# AVCS[[".hzr_status"]] <- with(AVCS, {
#   if (any(is.na(DEAD) | DEAD < 0)) {
#     stop("This job has rows where the EVENT count (DEAD) is missing or ",
#          "negative. ...")   # message trimmed here; see below
#   }
#   ifelse(DEAD > 0, 1, 0)
# })
# 
# # chunk: fit
# fit <- hazard(
#   data = AVCS, time = INT_DEAD, status = .hzr_status, fit = TRUE,
#   dist = "multiphase",
#   phases = list(
#     hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1, fixed = "m"),
#     hzr_phase("constant")
#   ),
#   theta = c(
#     log(0.2361727), log(0.1512095), 1.438652, 1, log(0.0005436977)
#   ),
#   weights = ifelse(DEAD > 0, DEAD, 1),
#   control = list(conserve = TRUE)
# )

