## -----------------------------------------------------------------------------
#| include: false
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE)
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
# fit <- hazard(
#   ...,
#   control = list(
#     nocov      = TRUE,   # suppress covariance output
#     nocor      = TRUE,   # suppress correlation output
#     condition  = 14      # CONDITION= switch
#   )
# )


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
#   theta   = c(MUE = 0.3504743, THALF = 0.1905077, NU = 1.437416,
#               M = 1,           MUC   = 4.391673e-07),
#   control = list(
#     fix = c("M")   # FIXM → freeze M during optimization
#   ),
#   ...
# )


## -----------------------------------------------------------------------------
#| eval: false
# # Build design matrix from the AVC data set
# X <- data.matrix(avcs[, c("AGE", "COM_IV", "MAL", "OPMOS", "OP_AGE",
#                            "STATUS", "INC_SURG", "ORIFICE")])
# 
# # Starting values from SAS EARLY + CONSTANT blocks combined
# theta_start <- c(
#   AGE      = -0.03205774,
#   COM_IV   =  1.336675,
#   MAL      =  0.6872028,
#   OPMOS    = -0.01963377,
#   OP_AGE   =  0.0002086689,
#   STATUS   =  0.5169533,   # EARLY phase coefficient
#   INC_SURG =  1.375285,
#   ORIFICE  =  3.11765
# )
# 
# fit <- hazard(
#   time   = avcs$INT_DEAD,
#   status = avcs$DEAD,
#   x      = X,
#   theta  = theta_start,
#   dist   = "weibull"
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
#| eval: false
# # Assumed: avcs is a data.frame read from the AVC flat file
# # (same variables as the SAS DATA step)
# 
# avcs <- avcs |>
#   transform(
#     LN_AGE   = log(AGE),
#     LN_OPMOS = log(OPMOS),
#     LN_INC   = ifelse(is.na(INC_SURG), NA, log(INC_SURG + 1)),
#     LN_NYHA  = log(STATUS)
#   )
# 
# # Replace missing INC_SURG with column mean (mirrors PROC STANDARD REPLACE)
# avcs$INC_SURG[is.na(avcs$INC_SURG)] <- mean(avcs$INC_SURG, na.rm = TRUE)
# 
# X <- data.matrix(avcs[, c("AGE", "COM_IV", "MAL", "OPMOS", "OP_AGE",
#                            "STATUS", "INC_SURG", "ORIFICE")])
# 
# fit <- hazard(
#   time    = avcs$INT_DEAD,
#   status  = avcs$DEAD,
#   x       = X,
#   theta   = c(
#     # Hazard shape parameters (from PARMS)
#     MUE   = 0.3504743,
#     THALF = 0.1905077,
#     NU    = 1.437416,
#     M     = 1,
#     MUC   = 4.391673e-07,
#     # Covariate coefficients (from EARLY + CONSTANT blocks)
#     AGE      = -0.03205774,
#     COM_IV   =  1.336675,
#     MAL      =  0.6872028,
#     OPMOS    = -0.01963377,
#     OP_AGE   =  0.0002086689,
#     STATUS   =  0.5169533,
#     INC_SURG =  1.375285,
#     ORIFICE  =  3.11765
#   ),
#   dist    = "weibull",
#   control = list(
#     condition = 14,
#     quasi     = TRUE,
#     conserve  = TRUE,
#     fix       = c("M")   # FIXM
#   )
# )
# 
# fit


## -----------------------------------------------------------------------------
#| eval: false
#| label: predict-example
# # Linear predictor (X %*% theta)
# eta <- predict(fit, type = "linear_predictor")
# 
# # Hazard scale (exp(eta))
# hz  <- predict(fit, type = "hazard")

