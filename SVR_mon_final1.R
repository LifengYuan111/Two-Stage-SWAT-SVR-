# ---------------------------------------------------------
# SWAT+–SVR Hybrid Monthly Streamflow Model
# ---------------------------------------------------------
# Study Area:
#   Roubidoux Creek Watershed, Missouri, USA
#
# Model Structure:
#   Stage 1: Baseflow-proxy SVR
#   Stage 2: Streamflow SVR
#
# Validation Strategy:
#   Strict temporal validation
#   Calibration period = first 70% of simulation period
#   Validation period  = remaining 30%
#
# Inputs:
#   - SWAT+ monthly streamflow
#   - Monthly precipitation
#   - Antecedent precipitation index (API)
#   - Rolling 3-month precipitation (RR3)
#   - Drainage area
#   - Seasonal sinusoidal predictors
#
# Output:
#   - SWAT+ baseline performance
#   - SWAT+SVR corrected streamflow
#   - Site-level and watershed-level metrics
#
# Author:
#   Lifeng Yuan
# ---------------------------------------------------------

# -----------------------------
# 0. Packages
# -----------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(scales)
library(e1071)
library(hydroGOF)
library(FlowScreen)
library(tibble)
library(patchwork)

# -----------------------------
# 1. User settings
# -----------------------------
set.seed(123)

min_valid_days_month <- 10
bf_a <- 0.98
bf_BFImax <- 0.80

# API decay parameter for monthly precipitation
api_k <- 0.9

station_info <- tibble(
  site_no = c("06928420", "06928380", "06928330", "06928320", "06928300"),
  swat_unit = c("1", "30", "148", "134", "95"),
  area_km2 = c(706.6, 21.7, 13.4, 25.2, 427.6),
  station_name = c(
    "Roubidoux Creek at Polla Rd bl Ft. Leonard Wood",
    "Upper Smith Branch bl Eng. Ponds at FLW, MO",
    "Hurd Hollow Trib at FLW (Outfall 14)",
    "Musgrave Hollow Trib at FLW",
    "Roubidoux Creek above Fort Leonard Wood, MO"
  )
)

site_order <- c("06928420", "06928380", "06928330", "06928320", "06928300")

# -----------------------------
# 2. Helper functions
# -----------------------------
safe_scale_train_test <- function(train_df, test_df, cols) {
  centers <- sapply(train_df[, cols, drop = FALSE], mean, na.rm = TRUE)
  sds <- sapply(train_df[, cols, drop = FALSE], sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- 1
  
  train_scaled <- train_df
  test_scaled  <- test_df
  
  train_scaled[, cols] <- scale(train_df[, cols, drop = FALSE], center = centers, scale = sds)
  test_scaled[, cols]  <- scale(test_df[, cols, drop = FALSE], center = centers, scale = sds)
  
  list(
    train = train_scaled,
    test = test_scaled,
    center = centers,
    scale = sds
  )
}

calc_metrics_one <- function(obs, sim) {
  obs <- as.numeric(obs)
  sim <- as.numeric(sim)
  
  keep <- is.finite(obs) & is.finite(sim)
  obs <- obs[keep]
  sim <- sim[keep]
  
  if (length(obs) < 2) {
    return(list(
      NSE = NA_real_,
      KGE = NA_real_,
      KGElf = NA_real_,
      PBIAS = NA_real_,
      RMSE = NA_real_,
      R2 = NA_real_,
      n = length(obs)
    ))
  }
  
  list(
    NSE = hydroGOF::NSE(sim = sim, obs = obs),
    KGE = hydroGOF::KGE(sim = sim, obs = obs),
    KGElf = hydroGOF::KGElf(sim = sim, obs = obs),
    PBIAS = hydroGOF::pbias(sim = sim, obs = obs),
    RMSE = hydroGOF::rmse(sim = sim, obs = obs),
    R2 = cor(obs, sim)^2,
    n = length(obs)
  )
}

calc_metrics_tbl <- function(df, obs_col, sim_col) {
  m <- calc_metrics_one(df[[obs_col]], df[[sim_col]])
  tibble(
    NSE = m$NSE,
    KGE = m$KGE,
    KGElf = m$KGElf,
    PBIAS = m$PBIAS,
    RMSE = m$RMSE,
    R2 = m$R2,
    n = m$n
  )
}

do_one_site_eckhardt_obs <- function(df_site, flow_col = "obs_q_daily", a = 0.98, BFImax = 0.8) {
  Q_raw <- df_site[[flow_col]]
  Q_fill <- Q_raw
  
  if (any(is.na(Q_fill))) {
    idx <- which(!is.na(Q_fill))
    if (length(idx) >= 2) {
      Q_fill <- approx(idx, Q_fill[idx], xout = seq_along(Q_fill), rule = 2)$y
    }
  }
  
  if (all(is.na(Q_fill)) || length(na.omit(Q_fill)) < 2) {
    df_site$obs_baseflow_daily <- NA_real_
    return(df_site)
  }
  
  bf <- FlowScreen::bf_eckhardt(discharge = Q_fill, a = a, BFI = BFImax)
  bf <- as.numeric(bf)
  bf[is.na(Q_raw)] <- NA_real_
  
  df_site$obs_baseflow_daily <- bf
  df_site
}

choose_date_scale <- function(dates) {
  span_years <- as.numeric(max(dates, na.rm = TRUE) - min(dates, na.rm = TRUE)) / 365.25
  
  if (span_years <= 2) {
    list(breaks = "3 months", labels = "%Y-%m")
  } else if (span_years <= 4) {
    list(breaks = "6 months", labels = "%Y-%m")
  } else if (span_years <= 10) {
    list(breaks = "1 year", labels = "%Y")
  } else {
    list(breaks = "2 years", labels = "%Y")
  }
}

calc_api <- function(p, k = 0.9) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  
  if (length(p) == 0) return(out)
  if (all(is.na(p))) return(out)
  
  for (i in seq_along(p)) {
    p_i <- ifelse(is.na(p[i]), 0, p[i])
    if (i == 1) {
      out[i] <- p_i
    } else {
      prev <- ifelse(is.na(out[i - 1]), 0, out[i - 1])
      out[i] <- p_i + k * prev
    }
  }
  out
}

calc_rr3 <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  if (length(p) == 0) return(out)
  
  for (i in seq_along(p)) {
    idx <- max(1, i - 2):i
    out[i] <- sum(ifelse(is.na(p[idx]), 0, p[idx]))
  }
  out
}

# -----------------------------
# 3. Clean inputs and read SWAT+ monthly output
# -----------------------------

obs_split_clean <- obs_split %>%
  ungroup() %>%
  mutate(
    site_no = sprintf("%08d", as.integer(gsub("[^0-9]", "", as.character(site_no)))),
    Date = as.Date(Date),
    Discharge_cms = as.numeric(Discharge_cms),
    Month = floor_date(Date, unit = "month")
  )

pcp_monthly <- pcp_data_mapped %>%
  mutate(
    site_no = sprintf("%08d", as.integer(gsub("[^0-9]", "", as.character(site_no)))),
    Date = as.Date(Date),
    Precip_mm = as.numeric(Precip_mm),
    Month = floor_date(Date, unit = "month")
  ) %>%
  group_by(site_no, Month) %>%
  summarise(
    precip_mm = sum(Precip_mm, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Read SWAT+ monthly channel output
# -----------------------------

channel_flow_mon <- read.table(
  "C:/SWAT Projects/FLW_0107_backup/Scenarios/Default/TxtInOut/channel_sd_mon.txt",
  skip = 3,
  header = TRUE
)

col_names <- c(
  "jday", "mon", "day", "yr", "unit", "gis_id", "name",
  "area", "precip", "evap", "seep", "flo_stor", "sed_stor",
  "orgn_stor", "sedp_stor", "no3_stor", "solp_stor", "chla_stor",
  "nh3_stor", "no2_stor", "cbod_stor", "dox_stor", "san_stor",
  "sil_stor", "cla_stor", "sag_stor", "lag_stor", "grv_stor",
  "null1", "flo_in", "sed_in", "orgn_in", "sedp_in", "no3_in",
  "solp_in", "chla_in", "nh3_in", "no2_in", "cbod_in", "dox_in",
  "san_in", "sil_in", "cla_in", "sag_in", "lag_in", "grv_in",
  "null2", "flo_out", "sed_out", "orgn_out", "sedp_out", "no3_out",
  "solp_out", "chla_out", "nh3_out", "no2_out", "cbod_out", "dox_out",
  "san_out", "sil_out", "cla_out", "sag_out", "lag_out", "grv_out",
  "null3", "water_temp"
)

names(channel_flow_mon) <- col_names

channel_flow_mon <- channel_flow_mon %>%
  mutate(
    unit = as.character(unit),
    Month = as.Date(sprintf("%04d-%02d-01", yr, mon)),
    sim_q = as.numeric(flo_out)
  ) %>%
  select(unit, Month, sim_q)

# -----------------------------
# 4. Aggregate observed daily flow and baseflow to monthly
# -----------------------------

obs_daily_bf <- obs_split_clean %>%
  select(site_no, Date, Month, period, obs_q_daily = Discharge_cms) %>%
  arrange(site_no, Date) %>%
  group_by(site_no) %>%
  group_modify(
    ~ do_one_site_eckhardt_obs(
      .x,
      flow_col = "obs_q_daily",
      a = bf_a,
      BFImax = bf_BFImax
    )
  ) %>%
  ungroup()

obs_monthly <- obs_daily_bf %>%
  group_by(site_no, period, Month) %>%
  summarise(
    n_obs_valid = sum(!is.na(obs_q_daily)),
    n_bf_valid = sum(!is.na(obs_baseflow_daily)),
    obs_q = ifelse(
      n_obs_valid >= min_valid_days_month,
      mean(obs_q_daily, na.rm = TRUE),
      NA_real_
    ),
    obs_baseflow_q = ifelse(
      n_bf_valid >= min_valid_days_month,
      mean(obs_baseflow_daily, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

# -----------------------------
# 5. Build monthly modeling dataset
# -----------------------------

sim_monthly <- station_info %>%
  mutate(
    site_no = sprintf("%08d", as.integer(gsub("[^0-9]", "", as.character(site_no)))),
    swat_unit = as.character(swat_unit)
  ) %>%
  select(site_no, station_name, swat_unit, area_km2) %>%
  left_join(
    channel_flow_mon,
    by = c("swat_unit" = "unit")
  )

monthly_all <- sim_monthly %>%
  left_join(
    obs_monthly,
    by = c("site_no", "Month")
  ) %>%
  left_join(
    pcp_monthly,
    by = c("site_no", "Month")
  ) %>%
  mutate(
    precip_mm = ifelse(is.na(precip_mm), 0, precip_mm)
  ) %>%
  filter(period %in% c("calibration", "validation")) %>%
  arrange(site_no, Month)


# -----------------------------
# 6. Derive hydrologic memory variables
# -----------------------------
# API represents antecedent wetness conditions using an
# exponentially decaying precipitation memory.
#
# RR3 represents cumulative precipitation over the current
# and previous two months.
#
# Lagged precipitation, lagged SWAT+ flow, lagged API, and
# lagged RR3 are included to represent delayed watershed
# response and hydrologic memory, which are important in
# karst systems.
# -----------------------------
monthly_all <- monthly_all %>%
  group_by(site_no) %>%
  arrange(Month, .by_group = TRUE) %>%
  mutate(
    api = calc_api(precip_mm, k = api_k),
    rr3 = calc_rr3(precip_mm),
    precip_lag1 = lag(precip_mm, 1),
    precip_lag2 = lag(precip_mm, 2),
    sim_lag1 = lag(sim_q, 1),
    sim_lag2 = lag(sim_q, 2),
    api_lag1 = lag(api, 1),
    rr3_lag1 = lag(rr3, 1),
    month_num = month(Month),
    sin12 = sin(2 * pi * month_num / 12),
    cos12 = cos(2 * pi * month_num / 12)
  ) %>%
  ungroup()

# -----------------------------
# 7. Split calibration / validation
# Strict validation:
# - observed baseflow used only as target in calibration
# - validation never uses observed baseflow directly as predictor
# -----------------------------
bf_train <- monthly_all %>%
  filter(period == "calibration") %>%
  filter(
    !is.na(obs_baseflow_q),
    !is.na(sim_q),
    !is.na(area_km2),
    !is.na(precip_mm),
    !is.na(precip_lag1),
    !is.na(precip_lag2),
    !is.na(api),
    !is.na(api_lag1),
    !is.na(rr3),
    !is.na(rr3_lag1),
    !is.na(sim_lag1),
    !is.na(sim_lag2)
  )

bf_all_pred <- monthly_all %>%
  filter(
    !is.na(sim_q),
    !is.na(area_km2),
    !is.na(precip_mm),
    !is.na(precip_lag1),
    !is.na(precip_lag2),
    !is.na(api),
    !is.na(api_lag1),
    !is.na(rr3),
    !is.na(rr3_lag1),
    !is.na(sim_lag1),
    !is.na(sim_lag2)
  )

# -----------------------------
# 8. Stage 1: Baseflow-proxy SVR feature preparation
# -----------------------------
# The target variable is monthly observed baseflow estimated
# from observed streamflow using the Eckhardt digital filter.
#
# The predictors include SWAT+-simulated flow, precipitation,
# API, RR3, lagged hydrologic variables, drainage area, and
# seasonal terms.
#
# A log1p transformation is used to reduce skewness and avoid
# problems with zero flow values.
# -----------------------------
make_bf_features <- function(df) {
  df %>%
    mutate(
      bf_y = log1p(pmax(obs_baseflow_q, 0)),
      bf_x_sim = log1p(pmax(sim_q, 0)),
      bf_x_p = log1p(pmax(precip_mm, 0)),
      bf_x_p1 = log1p(pmax(precip_lag1, 0)),
      bf_x_p2 = log1p(pmax(precip_lag2, 0)),
      bf_x_api = log1p(pmax(api, 0)),
      bf_x_api1 = log1p(pmax(api_lag1, 0)),
      bf_x_rr3 = log1p(pmax(rr3, 0)),
      bf_x_rr31 = log1p(pmax(rr3_lag1, 0)),
      bf_x_s1 = log1p(pmax(sim_lag1, 0)),
      bf_x_s2 = log1p(pmax(sim_lag2, 0)),
      bf_x_a = log1p(pmax(area_km2, 0))
    )
}

bf_train <- make_bf_features(bf_train)
bf_all_pred <- make_bf_features(bf_all_pred)

bf_feature_cols <- c(
  "bf_x_sim", "bf_x_p", "bf_x_p1", "bf_x_p2",
  "bf_x_api", "bf_x_api1",
  "bf_x_rr3", "bf_x_rr31",
  "bf_x_s1", "bf_x_s2", "bf_x_a",
  "sin12", "cos12"
)

if (nrow(bf_train) < 10) {
  stop("Too few rows in bf_train to fit baseflow proxy model.")
}

bf_scaled <- safe_scale_train_test(
  train_df = bf_train,
  test_df = bf_all_pred,
  cols = bf_feature_cols
)

bf_train_scaled <- bf_scaled$train
bf_all_scaled <- bf_scaled$test

# -----------------------------
# 9. Tune and fit baseflow proxy SVR
# -----------------------------
# set.seed(123)
# Option A: use best model from grid search
# bf_tune <- tune(
#   svm,
#   bf_y ~ bf_x_sim + bf_x_p + bf_x_p1 + bf_x_p2 +
#     bf_x_api + bf_x_api1 +
#     bf_x_rr3 + bf_x_rr31 +
#     bf_x_s1 + bf_x_s2 + bf_x_a +
#     sin12 + cos12,
#   data = as.data.frame(bf_train_scaled),
#   type = "eps-regression",
#   kernel = "radial",
#   ranges = list(
#     epsilon = 2^seq(-8, -1, by = 1),
#     gamma   = 2^seq(-8,  2, by = 1),
#     cost    = 2^seq(-1,  8, by = 1)
#   ),
#   tunecontrol = tune.control(sampling = "cross", cross = 5),
#   scale = FALSE
# )
# 
# bf_model <- bf_tune$best.model
# 
# if (is.null(bf_model)) {
#   stop("bf_model is NULL. Baseflow proxy model fitting failed.")
# }

bf_formula <-   bf_y ~ bf_x_sim + bf_x_p + bf_x_p1 + bf_x_p2 +
  bf_x_api + bf_x_api1 +
  bf_x_rr3 + bf_x_rr31 +
  bf_x_s1 + bf_x_s2 + bf_x_a +
  sin12 + cos12
# Option B: use fixed parameters selected from previous tuning
bf_model <- svm(
  formula = bf_formula,
  data = as.data.frame(bf_train_scaled),
  
  type = "eps-regression",
  kernel = "radial",
  
  epsilon = 0.125,
  gamma = 0.0625,
  cost = 16,
  
  scale = FALSE
)

bf_newdata <- bf_all_scaled %>%
  select(all_of(bf_feature_cols)) %>%
  as.data.frame()

bf_pred <- predict(bf_model, newdata = bf_newdata)
bf_pred <- as.numeric(bf_pred)

if (length(bf_pred) != nrow(bf_all_scaled)) {
  stop(
    paste0(
      "Baseflow proxy prediction length mismatch: predicted ",
      length(bf_pred),
      " rows, but bf_all_scaled has ",
      nrow(bf_all_scaled),
      " rows."
    )
  )
}

bf_all_scaled <- bf_all_scaled %>%
  mutate(
    bf_proxy_q = pmax(expm1(bf_pred), 0)
  )

monthly_model <- monthly_all %>%
  left_join(
    bf_all_scaled %>%
      select(site_no, Month, bf_proxy_q),
    by = c("site_no", "Month")
  )


# -----------------------------
# 10. Stage 2: Streamflow SVR feature preparation
# -----------------------------
# The final streamflow SVR predicts observed monthly streamflow
# using only two predictors: SWAT+-simulated monthly flow and
# the Stage-1 predicted baseflow proxy.
#
# In strict validation, observed baseflow is not used directly.
# The validation-period baseflow input is generated by the
# calibrated Stage 1 SVR model.
# -----------------------------

flow_data <- monthly_model %>%
  filter(period %in% c("calibration", "validation")) %>%
  filter(
    !is.na(obs_q),
    !is.na(sim_q),
    !is.na(bf_proxy_q)
  )

flow_train <- flow_data %>%
  filter(period == "calibration")

flow_test <- flow_data %>%
  filter(period == "validation")

make_flow_features_simple <- function(df) {
  df %>%
    mutate(
      y = log1p(pmax(obs_q, 0)),
      x_sim = log1p(pmax(sim_q, 0)),
      x_bf = log1p(pmax(bf_proxy_q, 0))
    )
}

flow_train <- make_flow_features_simple(flow_train)
flow_test  <- make_flow_features_simple(flow_test)

flow_feature_cols <- c("x_sim", "x_bf")

flow_scaled <- safe_scale_train_test(
  train_df = flow_train,
  test_df = flow_test,
  cols = flow_feature_cols
)

train_scaled <- flow_scaled$train
test_scaled  <- flow_scaled$test

# ---------------------------------------------------------
# 11. Tune and fit simple flow SVR
# ---------------------------------------------------------

set.seed(123)

flow_formula <- y ~ x_sim + x_bf

# flow_tune <- tune(
#   svm,
#   flow_formula,
#   data = as.data.frame(train_scaled),
#   type = "eps-regression",
#   kernel = "radial",
#   ranges = list(
#     epsilon = 2^seq(-8, -1, by = 1),
#     gamma   = 2^seq(-8,  2, by = 1),
#     cost    = 2^seq(-1,  8, by = 1)
#   ),
#   tunecontrol = tune.control(sampling = "cross", cross = 5),
#   scale = FALSE
# )
# 
# flow_model <- flow_tune$best.model
# 
# if (is.null(flow_model)) {
#   stop("flow_model is NULL. Simple flow SVR fitting failed.")
# }

flow_model <- svm(
  formula = flow_formula,
  data = as.data.frame(train_scaled),
  
  type = "eps-regression",
  kernel = "radial",
  
  epsilon = 0.00390626,
  gamma = 0.5,
  cost = 32,
  
  scale = FALSE
)


# ---------------------------------------------------------
# 12. Predict
# ---------------------------------------------------------

flow_new_train <- train_scaled %>%
  select(all_of(flow_feature_cols)) %>%
  as.data.frame()

flow_new_test <- test_scaled %>%
  select(all_of(flow_feature_cols)) %>%
  as.data.frame()

train_scaled$pred_svr <- expm1(predict(flow_model, newdata = flow_new_train))
test_scaled$pred_svr  <- expm1(predict(flow_model, newdata = flow_new_test))

train_scaled$pred_svr <- pmax(train_scaled$pred_svr, 0)
test_scaled$pred_svr  <- pmax(test_scaled$pred_svr, 0)

train_scaled$pred_swat <- train_scaled$sim_q
test_scaled$pred_swat  <- test_scaled$sim_q

# -----------------------------
# 13. Model performance evaluation
# -----------------------------
# Performance is evaluated separately for calibration and
# validation periods.
#
# Metrics include:
#   NSE    : overall predictive efficiency
#   KGE    : combined correlation, bias, and variability skill
#   KGElf  : low-flow focused KGE
#   PBIAS  : mean bias direction and magnitude
#   RMSE   : absolute simulation error
#   R2     : linear association between observed and simulated flow
#
# Site-level metrics are calculated independently for each gage.
# Overall metrics are calculated using pooled observations across
# all stations.
# -----------------------------

metrics_overall_simple <- bind_rows(
  calc_metrics_tbl(train_scaled, "obs_q", "pred_swat") %>%
    mutate(model = "SWAT+", period = "calibration"),
  
  calc_metrics_tbl(train_scaled, "obs_q", "pred_svr") %>%
    mutate(model = "SWAT+SVR_simple", period = "calibration"),
  
  calc_metrics_tbl(test_scaled, "obs_q", "pred_swat") %>%
    mutate(model = "SWAT+", period = "validation"),
  
  calc_metrics_tbl(test_scaled, "obs_q", "pred_svr") %>%
    mutate(model = "SWAT+SVR_simple", period = "validation")
) %>%
  select(period, model, NSE, KGE, KGElf, PBIAS, RMSE, R2, n)

print(metrics_overall_simple)

metrics_by_site_simple <- bind_rows(
  train_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_tbl(.x, "obs_q", "pred_swat")) %>%
    ungroup() %>%
    mutate(model = "SWAT+", period = "calibration"),
  
  train_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_tbl(.x, "obs_q", "pred_svr")) %>%
    ungroup() %>%
    mutate(model = "SWAT+SVR_simple", period = "calibration"),
  
  test_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_tbl(.x, "obs_q", "pred_swat")) %>%
    ungroup() %>%
    mutate(model = "SWAT+", period = "validation"),
  
  test_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_tbl(.x, "obs_q", "pred_svr")) %>%
    ungroup() %>%
    mutate(model = "SWAT+SVR_simple", period = "validation")
) %>%
  select(site_no, station_name, period, model, NSE, KGE, KGElf, PBIAS, RMSE, R2, n)

print(metrics_by_site_simple, n = 30)

# ---------------------------------------------------------
# 14. Prepare plotting data
# ---------------------------------------------------------

plot_base_simple <- bind_rows(train_scaled, test_scaled) %>%
  select(
    site_no, station_name, Month, period,
    obs_q, sim_q, pred_swat, pred_svr, precip_mm, bf_proxy_q
  ) %>%
  mutate(
    pred_swat = sim_q
  )

plot_df_simple <- plot_base_simple %>%
  select(site_no, station_name, Month, period, obs_q, pred_swat, pred_svr) %>%
  pivot_longer(
    cols = c(obs_q, pred_swat, pred_svr),
    names_to = "series",
    values_to = "flow_cms"
  ) %>%
  mutate(
    series = recode(
      series,
      obs_q = "Observed",
      pred_swat = "SWAT+",
      pred_svr = "SWAT+SVR_simple"
    )
  )

# -----------------------------
# 15. Plot observed, SWAT+, and SWAT+SVR monthly flow
# -----------------------------
# Each panel includes:
#   - monthly precipitation as an inverted bar plot
#   - observed monthly streamflow
#   - baseline SWAT+ monthly streamflow
#   - corrected SWAT+SVR monthly streamflow
#   - performance metrics for the selected period
# -----------------------------

plot_monthly_hybrid_simple <- function(site_id, target_period = c("calibration", "validation")) {
  
  target_period <- match.arg(target_period)
  
  df_flow <- plot_df_simple %>%
    filter(site_no == site_id, period == target_period) %>%
    arrange(Month)
  
  df_p <- plot_base_simple %>%
    select(site_no, station_name, Month, period, precip_mm) %>%
    distinct() %>%
    filter(site_no == site_id, period == target_period) %>%
    arrange(Month)
  
  if (nrow(df_flow) == 0) {
    return(
      ggplot() +
        annotate(
          "text", x = 1, y = 1,
          label = paste("No", target_period, "data for", site_id),
          size = 5
        ) +
        theme_void()
    )
  }
  
  station_lab <- unique(df_flow$station_name)[1]
  
  df_metric <- df_flow %>%
    select(Month, series, flow_cms) %>%
    group_by(Month, series) %>%
    summarise(
      flow_cms = mean(flow_cms, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = series, values_from = flow_cms) %>%
    arrange(Month)
  
  m_swat <- calc_metrics_one(
    obs = df_metric$Observed,
    sim = df_metric$`SWAT+`
  )
  
  m_svr <- calc_metrics_one(
    obs = df_metric$Observed,
    sim = df_metric$`SWAT+SVR_simple`
  )
  
  metrics_text <- paste0(
    "SWAT+             : ",
    "NSE=", sprintf("%.2f", m_swat$NSE),
    " | KGE=", sprintf("%.2f", m_swat$KGE),
    " | KGElf=", sprintf("%.2f", m_swat$KGElf),
    " | PBIAS=", sprintf("%.2f", m_swat$PBIAS), "%",
    " | RMSE=", sprintf("%.2f", m_swat$RMSE),
    " | n=", m_swat$n,
    "\n",
    "SWAT+SVR_simple   : ",
    "NSE=", sprintf("%.2f", m_svr$NSE),
    " | KGE=", sprintf("%.2f", m_svr$KGE),
    " | KGElf=", sprintf("%.2f", m_svr$KGElf),
    " | PBIAS=", sprintf("%.2f", m_svr$PBIAS), "%",
    " | RMSE=", sprintf("%.2f", m_svr$RMSE),
    " | n=", m_svr$n
  )
  
  max_p <- max(df_p$precip_mm, na.rm = TRUE)
  if (!is.finite(max_p) || max_p <= 0) max_p <- 1
  
  x_scale_i <- choose_date_scale(df_flow$Month)
  
  p1 <- ggplot(df_p, aes(x = Month, y = precip_mm)) +
    geom_col(fill = "#A6CEE3", width = 25) +
    scale_y_reverse(
      limits = c(max_p * 1.05, 0),
      expand = c(0, 0),
      sec.axis = sec_axis(~ ., name = "Monthly precipitation (mm)")
    ) +
    scale_x_date(
      date_breaks = x_scale_i$breaks,
      labels = date_format(x_scale_i$labels),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.y.left = element_blank(),
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.title.y.right = element_text(size = 11),
      axis.text.y.right = element_text(size = 9),
      plot.margin = margin(2, 100, 0, 20)
    )
  
  p2 <- ggplot(df_flow, aes(x = Month, y = flow_cms, color = series)) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Observed" = "#D94801",
        "SWAT+" = "#1F78B4",
        "SWAT+SVR_simple" = "#33A02C"
      ),
      name = NULL
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = metrics_text,
      hjust = 1.02,
      vjust = 1.05,
      size = 3.0,
      lineheight = 1.05
    ) +
    scale_x_date(
      date_breaks = x_scale_i$breaks,
      labels = date_format(x_scale_i$labels),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = paste0(site_id, " | ", station_lab, " | ", target_period),
      x = "Month",
      y = expression(paste("Flow (m"^3, " s"^-1, ")"))
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      legend.justification = "center",
      plot.title = element_text(face = "bold", size = 11),
      plot.margin = margin(0, 100, 8, 20)
    )
  
  p1 / p2
}

# ---------------------------------------------------------
# 16. Full plots
# ---------------------------------------------------------

plot_list_cal_simple <- lapply(site_order, function(sid) {
  plot_monthly_hybrid_simple(sid, "calibration")
})

p_cal_full_simple <- wrap_plots(plot_list_cal_simple, ncol = 1, guides = "collect") &
  theme(
    legend.position = "top",
    legend.justification = "center"
  )

plot_list_val_simple <- lapply(site_order, function(sid) {
  plot_monthly_hybrid_simple(sid, "validation")
})

p_val_full_simple <- wrap_plots(plot_list_val_simple, ncol = 1, guides = "collect") &
  theme(
    legend.position = "top",
    legend.justification = "center"
  )

print(p_cal_full_simple)
print(p_val_full_simple)

# ---------------------------------------------------------
# 17. Optional single-site plots
# ---------------------------------------------------------

p_cal_28420_simple <- plot_monthly_hybrid_simple("06928420", "calibration")
p_val_28420_simple <- plot_monthly_hybrid_simple("06928420", "validation")

print(p_cal_28420_simple)
print(p_val_28420_simple)

# -----------------------------
# 18. Export outputs
# -----------------------------
# The following files are exported:
#   1. Overall performance metrics
#   2. Site-level performance metrics
#   3. Monthly plotting dataset
#   4. Best SVR hyperparameters
#   5. Calibration and validation figures
# -----------------------------

output_dir_simple <- "C:/SWAT Projects/Calibration/Monthly_Hybrid/Simple_QSWAT_BFhat"

dir.create(output_dir_simple, recursive = TRUE, showWarnings = FALSE)

write.csv(
  metrics_overall_simple,
  file.path(output_dir_simple, "metrics_overall_simple_QSWAT_BFhat.csv"),
  row.names = FALSE
)

write.csv(
  metrics_by_site_simple,
  file.path(output_dir_simple, "metrics_by_site_simple_QSWAT_BFhat.csv"),
  row.names = FALSE
)

write.csv(
  plot_base_simple,
  file.path(output_dir_simple, "monthly_plot_data_simple_QSWAT_BFhat.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(flow_tune$best.parameters),
  file.path(output_dir_simple, "monthly_stage2_simple_flow_best_parameters.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(output_dir_simple, "monthly_simple_calibration.tiff"),
  plot = p_cal_full_simple,
  width = 12,
  height = 18,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir_simple, "monthly_simple_validation.tiff"),
  plot = p_val_full_simple,
  width = 12,
  height = 18,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir_simple, "monthly_simple_calibration.pdf"),
  plot = p_cal_full_simple,
  width = 12,
  height = 18,
  units = "in"
)

ggsave(
  filename = file.path(output_dir_simple, "monthly_simple_validation.pdf"),
  plot = p_val_full_simple,
  width = 12,
  height = 18,
  units = "in"
)

cat("\n====================================================\n")
cat("Simple monthly SWAT+SVR model finished.\n")
cat("Stage 2 predictors: Q_SWAT + BF_hat only\n")
cat("Outputs saved to:\n")
cat(output_dir_simple, "\n")
cat("====================================================\n")
