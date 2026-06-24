# =========================================================
# Supplementary Code S2
# Daily Two-Stage SWAT+–SVR Hybrid Streamflow Model
#
# Study area:
#   Roubidoux Creek Watershed, Missouri, USA
#
# Model structure:
#   Stage 1: SVR-based baseflow proxy estimation
#   Stage 2: SVR-based daily streamflow correction
#
# Validation strategy:
#   Observed baseflow is derived only during calibration.
#   Validation-period baseflow is predicted by the Stage 1 SVR.
#   This avoids information leakage from validation observations.
#
# Required input objects:
#   obs_split        : observed daily streamflow
#   channel_flow     : SWAT+ daily channel output
#   pcp_data_mapped  : daily precipitation mapped to each site
# =========================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(scales)
library(e1071)
library(hydroGOF)
library(FlowScreen)
library(tibble)
library(patchwork)

set.seed(123)

# ---------------------------------------------------------
# 1. User-defined settings
# ---------------------------------------------------------

# Eckhardt digital-filter parameters.
# a controls groundwater recession persistence.
# BFImax controls the maximum long-term baseflow fraction.
bf_a <- 0.98
bf_BFImax <- 0.80

# Antecedent Precipitation Index decay coefficient.
api_k <- 0.95

# Log-transform non-negative hydrologic variables.
use_log_transform_baseflow <- TRUE
use_log_transform_flow <- TRUE

output_dir <- "C:/SWAT Projects/Calibration/Daily_Hybrid_TwoStage"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Fixed SVR hyperparameters selected from prior tuning.
bf_svr_params <- list(
  epsilon = 0.125,
  gamma = 0.25,
  cost = 8
)

flow_svr_params <- list(
  epsilon = 0.25,
  gamma = 0.03125,
  cost = 32
)

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

site_order <- station_info$site_no

# ---------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------

safe_scale_train_test <- function(train_df, test_df, cols) {
  # Standardize predictors using calibration-period statistics only.
  # The same scaling parameters are applied to validation data.
  
  centers <- sapply(train_df[, cols, drop = FALSE], mean, na.rm = TRUE)
  sds <- sapply(train_df[, cols, drop = FALSE], sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- 1
  
  train_scaled <- train_df
  test_scaled <- test_df
  
  train_scaled[, cols] <- scale(train_df[, cols, drop = FALSE],
                                center = centers, scale = sds)
  test_scaled[, cols] <- scale(test_df[, cols, drop = FALSE],
                               center = centers, scale = sds)
  
  list(train = train_scaled, test = test_scaled)
}

calc_metrics_one <- function(obs, sim) {
  # Hydrologic model performance metrics.
  
  obs <- as.numeric(obs)
  sim <- as.numeric(sim)
  keep <- is.finite(obs) & is.finite(sim)
  obs <- obs[keep]
  sim <- sim[keep]
  
  if (length(obs) < 2) {
    return(tibble(
      NSE = NA_real_, KGE = NA_real_, KGElf = NA_real_,
      PBIAS = NA_real_, RMSE = NA_real_, R2 = NA_real_,
      n = length(obs)
    ))
  }
  
  tibble(
    NSE = hydroGOF::NSE(sim = sim, obs = obs),
    KGE = hydroGOF::KGE(sim = sim, obs = obs),
    KGElf = tryCatch(
      hydroGOF::KGElf(sim = sim, obs = obs),
      error = function(e) NA_real_
    ),
    PBIAS = hydroGOF::pbias(sim = sim, obs = obs),
    RMSE = hydroGOF::rmse(sim = sim, obs = obs),
    R2 = cor(obs, sim)^2,
    n = length(obs)
  )
}

do_one_site_eckhardt <- function(df_site, flow_col, out_col,
                                 a = 0.98, BFImax = 0.8) {
  # Apply Eckhardt digital filter to one site.
  # Missing values are temporarily interpolated for filter continuity
  # and then restored to NA.
  
  Q_raw <- df_site[[flow_col]]
  Q_fill <- Q_raw
  
  if (any(is.na(Q_fill))) {
    idx <- which(!is.na(Q_fill))
    if (length(idx) >= 2) {
      Q_fill <- approx(idx, Q_fill[idx],
                       xout = seq_along(Q_fill), rule = 2)$y
    }
  }
  
  if (all(is.na(Q_fill)) || length(na.omit(Q_fill)) < 2) {
    df_site[[out_col]] <- NA_real_
    return(df_site)
  }
  
  bf <- tryCatch(
    FlowScreen::bf_eckhardt(discharge = Q_fill, a = a, BFI = BFImax),
    error = function(e) rep(NA_real_, length(Q_fill))
  )
  
  bf <- as.numeric(bf)
  bf[is.na(Q_raw)] <- NA_real_
  df_site[[out_col]] <- bf
  
  df_site
}

calc_api <- function(p, k = 0.95) {
  # Antecedent Precipitation Index:
  # API_t = P_t + k * API_{t-1}
  
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  
  if (length(p) == 0) return(out)
  
  out[1] <- ifelse(is.na(p[1]), 0, p[1])
  
  if (length(p) >= 2) {
    for (i in 2:length(p)) {
      out[i] <- ifelse(is.na(p[i]), 0, p[i]) + k * out[i - 1]
    }
  }
  
  out
}

# ---------------------------------------------------------
# 3. Clean input data
# ---------------------------------------------------------

obs_split <- obs_split %>%
  ungroup() %>%
  mutate(
    site_no = sprintf("%08d", as.integer(gsub("[^0-9]", "", as.character(site_no)))),
    Date = as.Date(Date),
    Discharge_cms = as.numeric(Discharge_cms),
    period = as.character(period)
  )

channel_flow <- channel_flow %>%
  mutate(
    unit = as.character(unit),
    Date = as.Date(Date),
    flo_out = as.numeric(flo_out)
  )

pcp_data_mapped <- pcp_data_mapped %>%
  mutate(
    site_no = sprintf("%08d", as.integer(gsub("[^0-9]", "", as.character(site_no)))),
    Date = as.Date(Date),
    Precip_mm = as.numeric(Precip_mm)
  ) %>%
  select(site_no, Date, Precip_mm)

# ---------------------------------------------------------
# 4. Build daily modeling dataset
# ---------------------------------------------------------

sim_daily <- station_info %>%
  select(site_no, swat_unit, area_km2, station_name) %>%
  left_join(
    channel_flow %>% select(unit, Date, flo_out),
    by = c("swat_unit" = "unit")
  ) %>%
  rename(sim_q = flo_out)

obs_daily <- obs_split %>%
  select(site_no, Date, obs_q = Discharge_cms, period)

daily_all <- sim_daily %>%
  left_join(obs_daily, by = c("site_no", "Date")) %>%
  left_join(pcp_data_mapped, by = c("site_no", "Date")) %>%
  mutate(
    precip_mm = ifelse(is.na(Precip_mm), 0, Precip_mm)
  ) %>%
  select(site_no, swat_unit, area_km2, station_name,
         Date, sim_q, obs_q, period, precip_mm) %>%
  arrange(site_no, Date)

# ---------------------------------------------------------
# 5. Derive daily hydrologic predictors
# ---------------------------------------------------------

daily_all <- daily_all %>%
  group_by(site_no) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    api = calc_api(precip_mm, k = api_k),
    precip_lag1 = lag(precip_mm, 1),
    sim_q_lag1 = lag(sim_q, 1),
    doy = yday(Date),
    year_days = ifelse(leap_year(Date), 366, 365),
    sin_doy = sin(2 * pi * doy / year_days),
    cos_doy = cos(2 * pi * doy / year_days)
  ) %>%
  ungroup()

data_all <- daily_all %>%
  filter(period %in% c("calibration", "validation")) %>%
  filter(!is.na(obs_q)) %>%
  arrange(site_no, Date)

train_data <- data_all %>% filter(period == "calibration")
test_data <- data_all %>% filter(period == "validation")

if (nrow(train_data) == 0) stop("No calibration data found.")
if (nrow(test_data) == 0) stop("No validation data found.")

# ---------------------------------------------------------
# 6. Calibration-only observed baseflow separation
# ---------------------------------------------------------

train_data <- train_data %>%
  group_by(site_no) %>%
  group_modify(
    ~ do_one_site_eckhardt(
      .x,
      flow_col = "obs_q",
      out_col = "obs_baseflow_q",
      a = bf_a,
      BFImax = bf_BFImax
    )
  ) %>%
  ungroup()

# Validation-period observed baseflow is intentionally not used.
test_data$obs_baseflow_q <- NA_real_

# ---------------------------------------------------------
# 7. Stage 1: Baseflow-proxy SVR
# ---------------------------------------------------------

make_baseflow_features <- function(df) {
  df %>%
    mutate(
      y_bf = if (use_log_transform_baseflow)
        log1p(pmax(obs_baseflow_q, 0)) else obs_baseflow_q,
      
      bf_x_sim = if (use_log_transform_baseflow)
        log1p(pmax(sim_q, 0)) else sim_q,
      
      bf_x_p = if (use_log_transform_baseflow)
        log1p(pmax(precip_mm, 0)) else precip_mm,
      
      bf_x_p1 = if (use_log_transform_baseflow)
        log1p(pmax(precip_lag1, 0)) else precip_lag1,
      
      bf_x_api = if (use_log_transform_baseflow)
        log1p(pmax(api, 0)) else api,
      
      bf_x_area = if (use_log_transform_baseflow)
        log1p(pmax(area_km2, 0)) else area_km2
    )
}

train_bf <- make_baseflow_features(train_data)
test_bf <- make_baseflow_features(test_data)

baseflow_feature_cols <- c(
  "bf_x_sim",
  "bf_x_p",
  "bf_x_p1",
  "bf_x_api",
  "bf_x_area",
  "sin_doy",
  "cos_doy"
)

missing_vars <- setdiff(baseflow_feature_cols, colnames(train_bf))
if (length(missing_vars) > 0) {
  stop(paste("Missing baseflow features:",
             paste(missing_vars, collapse = ", ")))
}

train_bf <- train_bf %>%
  drop_na(all_of(c("y_bf", baseflow_feature_cols)))

test_bf <- test_bf %>%
  drop_na(all_of(baseflow_feature_cols))

bf_scaled <- safe_scale_train_test(
  train_df = train_bf,
  test_df = test_bf,
  cols = baseflow_feature_cols
)

train_bf_scaled <- bf_scaled$train
test_bf_scaled <- bf_scaled$test

bf_formula <- as.formula(
  paste("y_bf ~", paste(baseflow_feature_cols, collapse = " + "))
)

bf_model <- svm(
  formula = bf_formula,
  data = as.data.frame(train_bf_scaled),
  type = "eps-regression",
  kernel = "radial",
  epsilon = bf_svr_params$epsilon,
  gamma = bf_svr_params$gamma,
  cost = bf_svr_params$cost,
  scale = FALSE
)

pred_bf_train_raw <- predict(
  bf_model,
  newdata = as.data.frame(train_bf_scaled[, baseflow_feature_cols])
)

pred_bf_test_raw <- predict(
  bf_model,
  newdata = as.data.frame(test_bf_scaled[, baseflow_feature_cols])
)

if (use_log_transform_baseflow) {
  train_bf_scaled$bf_proxy_pred <- expm1(pred_bf_train_raw)
  test_bf_scaled$bf_proxy_pred <- expm1(pred_bf_test_raw)
} else {
  train_bf_scaled$bf_proxy_pred <- pred_bf_train_raw
  test_bf_scaled$bf_proxy_pred <- pred_bf_test_raw
}

train_bf_scaled$bf_proxy_pred <- pmax(train_bf_scaled$bf_proxy_pred, 0)
test_bf_scaled$bf_proxy_pred <- pmax(test_bf_scaled$bf_proxy_pred, 0)

metrics_baseflow_train <- calc_metrics_one(
  train_bf_scaled$obs_baseflow_q,
  train_bf_scaled$bf_proxy_pred
) %>%
  mutate(period = "calibration", model = "Baseflow proxy SVR")

print(metrics_baseflow_train)

# ---------------------------------------------------------
# 8. Merge predicted baseflow proxy into flow dataset
# ---------------------------------------------------------

train_flow0 <- train_data %>%
  left_join(
    train_bf_scaled %>% select(site_no, Date, bf_proxy_pred),
    by = c("site_no", "Date")
  )

test_flow0 <- test_data %>%
  left_join(
    test_bf_scaled %>% select(site_no, Date, bf_proxy_pred),
    by = c("site_no", "Date")
  )

# ---------------------------------------------------------
# 9. Stage 2: Daily streamflow SVR
# ---------------------------------------------------------

flow_combined <- bind_rows(train_flow0, test_flow0) %>%
  group_by(site_no) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    bf_proxy_lag1 = lag(bf_proxy_pred, 1)
  ) %>%
  ungroup()

train_flow <- flow_combined %>% filter(period == "calibration")
test_flow <- flow_combined %>% filter(period == "validation")

make_flow_features <- function(df) {
  df %>%
    mutate(
      y = if (use_log_transform_flow)
        log1p(pmax(obs_q, 0)) else obs_q,
      
      x_sim = if (use_log_transform_flow)
        log1p(pmax(sim_q, 0)) else sim_q,
      
      x_sim1 = if (use_log_transform_flow)
        log1p(pmax(sim_q_lag1, 0)) else sim_q_lag1,
      
      x_bf = if (use_log_transform_flow)
        log1p(pmax(bf_proxy_pred, 0)) else bf_proxy_pred,
      
      x_bf1 = if (use_log_transform_flow)
        log1p(pmax(bf_proxy_lag1, 0)) else bf_proxy_lag1,
      
      x_p = if (use_log_transform_flow)
        log1p(pmax(precip_mm, 0)) else precip_mm,
      
      x_p1 = if (use_log_transform_flow)
        log1p(pmax(precip_lag1, 0)) else precip_lag1,
      
      x_api = if (use_log_transform_flow)
        log1p(pmax(api, 0)) else api,
      
      x_area = if (use_log_transform_flow)
        log1p(pmax(area_km2, 0)) else area_km2
    )
}

train_flow <- make_flow_features(train_flow)
test_flow <- make_flow_features(test_flow)

flow_feature_cols <- c(
  "x_sim",
  "x_sim1",
  "x_bf",
  "x_bf1",
  "x_p",
  "x_p1",
  "x_api",
  "x_area",
  "sin_doy",
  "cos_doy"
)

missing_flow_vars <- setdiff(flow_feature_cols, colnames(train_flow))
if (length(missing_flow_vars) > 0) {
  stop(paste("Missing flow features:",
             paste(missing_flow_vars, collapse = ", ")))
}

train_flow <- train_flow %>%
  drop_na(all_of(c("y", flow_feature_cols)))

test_flow <- test_flow %>%
  drop_na(all_of(c("y", flow_feature_cols)))

flow_scaled <- safe_scale_train_test(
  train_df = train_flow,
  test_df = test_flow,
  cols = flow_feature_cols
)

train_flow_scaled <- flow_scaled$train
test_flow_scaled <- flow_scaled$test

flow_formula <- as.formula(
  paste("y ~", paste(flow_feature_cols, collapse = " + "))
)

flow_model <- svm(
  formula = flow_formula,
  data = as.data.frame(train_flow_scaled),
  type = "eps-regression",
  kernel = "radial",
  epsilon = flow_svr_params$epsilon,
  gamma = flow_svr_params$gamma,
  cost = flow_svr_params$cost,
  scale = FALSE
)

pred_train_raw <- predict(
  flow_model,
  newdata = as.data.frame(train_flow_scaled[, flow_feature_cols])
)

pred_test_raw <- predict(
  flow_model,
  newdata = as.data.frame(test_flow_scaled[, flow_feature_cols])
)

if (use_log_transform_flow) {
  train_flow_scaled$pred_svr <- expm1(pred_train_raw)
  test_flow_scaled$pred_svr <- expm1(pred_test_raw)
} else {
  train_flow_scaled$pred_svr <- pred_train_raw
  test_flow_scaled$pred_svr <- pred_test_raw
}

train_flow_scaled$pred_svr <- pmax(train_flow_scaled$pred_svr, 0)
test_flow_scaled$pred_svr <- pmax(test_flow_scaled$pred_svr, 0)

train_flow_scaled$pred_swat <- train_flow_scaled$sim_q
test_flow_scaled$pred_swat <- test_flow_scaled$sim_q

# ---------------------------------------------------------
# 10. Model performance evaluation
# ---------------------------------------------------------

metrics_overall <- bind_rows(
  calc_metrics_one(train_flow_scaled$obs_q, train_flow_scaled$pred_swat) %>%
    mutate(period = "calibration", model = "SWAT+"),
  
  calc_metrics_one(train_flow_scaled$obs_q, train_flow_scaled$pred_svr) %>%
    mutate(period = "calibration", model = "SWAT+SVR"),
  
  calc_metrics_one(test_flow_scaled$obs_q, test_flow_scaled$pred_swat) %>%
    mutate(period = "validation", model = "SWAT+"),
  
  calc_metrics_one(test_flow_scaled$obs_q, test_flow_scaled$pred_svr) %>%
    mutate(period = "validation", model = "SWAT+SVR")
)

metrics_by_site <- bind_rows(
  train_flow_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_one(.x$obs_q, .x$pred_swat)) %>%
    ungroup() %>%
    mutate(period = "calibration", model = "SWAT+"),
  
  train_flow_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_one(.x$obs_q, .x$pred_svr)) %>%
    ungroup() %>%
    mutate(period = "calibration", model = "SWAT+SVR"),
  
  test_flow_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_one(.x$obs_q, .x$pred_swat)) %>%
    ungroup() %>%
    mutate(period = "validation", model = "SWAT+"),
  
  test_flow_scaled %>%
    group_by(site_no, station_name) %>%
    group_modify(~ calc_metrics_one(.x$obs_q, .x$pred_svr)) %>%
    ungroup() %>%
    mutate(period = "validation", model = "SWAT+SVR")
)

print(metrics_overall)
print(metrics_by_site, n = 24)

# ---------------------------------------------------------
# 11. Prepare plotting and export datasets
# ---------------------------------------------------------

plot_base_daily <- bind_rows(train_flow_scaled, test_flow_scaled) %>%
  select(
    site_no, station_name, Date, period,
    obs_q, sim_q, pred_swat, pred_svr,
    precip_mm, bf_proxy_pred, api
  ) %>%
  mutate(pred_swat = sim_q)

write.csv(
  metrics_baseflow_train,
  file.path(output_dir, "daily_stage1_baseflow_metrics.csv"),
  row.names = FALSE
)

write.csv(
  metrics_overall,
  file.path(output_dir, "daily_stage2_flow_metrics_overall.csv"),
  row.names = FALSE
)

write.csv(
  metrics_by_site,
  file.path(output_dir, "daily_stage2_flow_metrics_by_site.csv"),
  row.names = FALSE
)

write.csv(
  plot_base_daily,
  file.path(output_dir, "daily_hybrid_plot_data.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    model = c("Stage 1 baseflow proxy SVR", "Stage 2 streamflow SVR"),
    epsilon = c(bf_svr_params$epsilon, flow_svr_params$epsilon),
    gamma = c(bf_svr_params$gamma, flow_svr_params$gamma),
    cost = c(bf_svr_params$cost, flow_svr_params$cost)
  ),
  file.path(output_dir, "daily_svr_fixed_hyperparameters.csv"),
  row.names = FALSE
)

cat("\n====================================================\n")
cat("Daily two-stage SWAT+–SVR model completed.\n")
cat("Outputs saved to:\n")
cat(output_dir, "\n")
cat("====================================================\n")