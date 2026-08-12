
# CLAIM RESERVING WITH CALENDAR-BASED ROLLING-ORIGIN VALIDATION
# ============================================================
#
# Models:
# A. Quasi-Poisson GLM
# B. Tweedie GLM
# C. LASSO-penalized Poisson regression (L1)
# Each final model is refitted using all 105 observed cells
# and then used to predict the 91 unobserved cells in the
# lower triangle.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cplm)
  library(glmnet)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# ------------------------------------------------------------
# 1. Cumulative paid-claims triangle
# ------------------------------------------------------------
cum_triangle <- matrix(c(
  69672343,92381070,92538892,92631912,92631912,92631912,92631912,92631912,92631912,92631912,92631912,92631912,92631912,92631912,
  85370339,107165385,107311953,107311953,107311953,107311953,107311953,107311953,107311953,107311953,107311953,107311953,107311953,NA,
  99765261,127078378,127419050,127419050,127419050,127419050,127419050,127419050,127419050,127419050,127419050,127419050,NA,NA,
  87220584,126166049,126596965,126596965,126612485,126612485,126612485,126612485,126612485,126612485,126612485,NA,NA,NA,
  75609735,93993020,94210309,94232449,94232449,94232449,94232449,94232449,94232449,94232449,NA,NA,NA,NA,
  113677021,152600580,153166804,153167204,153167204,153167204,153167204,153167204,153167204,NA,NA,NA,NA,NA,
  144455980,189810239,190567087,190767578,190767578,191642554,191867205,191867205,NA,NA,NA,NA,NA,NA,
  114860146,151695727,151824949,152008875,152012875,152049659,152049659,NA,NA,NA,NA,NA,NA,NA,
  145025581,173843004,174926631,175061051,175099551,175099551,NA,NA,NA,NA,NA,NA,NA,NA,
  128691752,163689327,163824793,163911628,163911628,NA,NA,NA,NA,NA,NA,NA,NA,NA,
  141608514,169466442,169530731,169537131,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
  154641559,185335558,185684038,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
  178455867,208492055,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
  174183278,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
), nrow = 14L, byrow = TRUE)

N <- nrow(cum_triangle)
stopifnot(N == 14L, ncol(cum_triangle) == N)

origin_labels <- c(
  "2018_S1", "2018_S2", "2019_S1", "2019_S2",
  "2020_S1", "2020_S2", "2021_S1", "2021_S2",
  "2022_S1", "2022_S2", "2023_S1", "2023_S2",
  "2024_S1", "2024_S2"
)
calendar_labels <- setNames(origin_labels, as.character(seq_len(N)))

# ------------------------------------------------------------
# 2. Convert to incremental claims and long format
# ------------------------------------------------------------
inc_triangle <- cum_triangle
for (i in seq_len(N)) {
  for (j in seq_len(N)) {
    if (is.na(cum_triangle[i, j])) {
      inc_triangle[i, j] <- NA_real_
    } else if (j == 1L) {
      inc_triangle[i, j] <- cum_triangle[i, j]
    } else {
      inc_triangle[i, j] <- cum_triangle[i, j] - cum_triangle[i, j - 1L]
    }
  }
}

wide_inc <- as.data.table(inc_triangle)
wide_inc[, acc := .I]
long_dat <- melt(
  wide_inc, id.vars = "acc", variable.name = "dev", value.name = "pmts"
)
long_dat[, dev := as.integer(sub("^V", "", dev))]
long_dat[, cal := acc + dev - 1L]
long_dat[, observed := !is.na(pmts)]
long_dat[, origin_period := origin_labels[acc]]
setorder(long_dat, acc, dev)

obs_dat <- long_dat[observed == TRUE, .(pmts, acc, dev, cal)]
future_dat <- long_dat[observed == FALSE, .(acc, dev, cal)]
stopifnot(nrow(obs_dat) == 105L, nrow(future_dat) == 91L)

outer_cutoffs <- 8:13
inner_start_cutoff <- 5L

# Exactly the same split is used for every model.
# dev<=cutoff prevents testing a development level absent from training.
calendar_split <- function(data, cutoff) {
  list(
    train = data[cal <= cutoff],
    test = data[cal == cutoff + 1L & dev > 1L & dev <= cutoff]
  )
}

calculate_metrics <- function(predictions) {
  predictions[, {
    sse <- sum((actual - predicted)^2)
    sst <- sum((actual - mean(actual))^2)
    .(
      N_test = .N,
      RMSE = sqrt(mean((actual - predicted)^2)),
      MAE = mean(abs(actual - predicted)),
      R_squared_predictive = 1 - sse / sst,
      Actual_total = sum(actual),
      Predicted_total = sum(predicted),
      Aggregate_relative_error =
        (sum(predicted) - sum(actual)) / sum(actual)
    )
  }]
}

calculate_metrics_by_period <- function(predictions) {
  predictions[, .(
    N_test = .N,
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAE = mean(abs(actual - predicted)),
    Actual_total = sum(actual),
    Predicted_total = sum(predicted),
    Aggregate_relative_error =
      (sum(predicted) - sum(actual)) / sum(actual)
  ), by = .(cutoff, train_through, test_period)]
}

latest_observed <- data.table(
  acc = seq_len(N),
  origin_period = origin_labels,
  observed_cumulative_paid = apply(
    cum_triangle, 1L, function(x) tail(na.omit(x), 1L)
  )
)

make_reserve_tables <- function(model_name, future_predictions) {
  fp <- cbind(copy(future_dat), predicted_incremental = future_predictions)
  fp[, `:=`(model = model_name, origin_period = origin_labels[acc])]
  setorder(fp, acc, dev)

  by_origin <- fp[, .(
    reserve = sum(predicted_incremental),
    future_cells = .N
  ), by = .(model, acc, origin_period)]

  full_origins <- data.table(
    model = model_name, acc = seq_len(N), origin_period = origin_labels
  )
  by_origin <- merge(
    full_origins, by_origin,
    by = c("model", "acc", "origin_period"), all.x = TRUE
  )
  by_origin[is.na(reserve), `:=`(reserve = 0, future_cells = 0L)]
  by_origin <- merge(
    by_origin, latest_observed,
    by = c("acc", "origin_period"), all.x = TRUE
  )
  by_origin[, ultimate_claims := observed_cumulative_paid + reserve]
  setorder(by_origin, acc)

  totals <- by_origin[, .(
    observed_cumulative_paid = sum(observed_cumulative_paid),
    reserve = sum(reserve),
    ultimate_claims = sum(ultimate_claims)
  ), by = model]

  list(cell_predictions = fp, by_origin = by_origin, totals = totals)
}
# Align factor levels for prediction
prepare_glm_data <- function(train, test) {
  tr <- copy(train)
  te <- copy(test)
  tr[, `:=`(acc_f = factor(acc), dev_f = factor(dev))]
  te[, `:=`(
    acc_f = factor(acc, levels = levels(tr$acc_f)),
    dev_f = factor(dev, levels = levels(tr$dev_f))
  )]
  if (anyNA(te$acc_f) || anyNA(te$dev_f)) {
    stop("Test data contain an origin/development level absent from training.")
  }
  list(train = tr, test = te)
}

# ============================================================
# A. QUASI-POISSON GLM: COMPLETE ANALYSIS
# ============================================================
cat("\n========================================\n")
cat("A. QUASI-POISSON GLM\n")
cat("========================================\n")

fit_predict_quasipoisson <- function(train, test) {
  dat <- prepare_glm_data(train, test)
  fit <- glm(
    pmts ~ acc_f + dev_f,
    family = quasipoisson(link = "log"),
    data = dat$train,
    control = glm.control(maxit = 200L, epsilon = 1e-8)
  )
  pred <- pmax(as.numeric(predict(fit, newdata = dat$test, type = "response")), 0)
  if (length(pred) != nrow(test) || anyNA(pred) || any(!is.finite(pred))) {
    stop("Quasi-Poisson produced invalid predictions.")
  }
  list(prediction = pred, fit = fit)
}

qp_roll_list <- vector("list", length(outer_cutoffs))
for (h in seq_along(outer_cutoffs)) {
  cutoff <- outer_cutoffs[h]
  split <- calendar_split(obs_dat, cutoff)
  ans <- fit_predict_quasipoisson(split$train, split$test)
  qp_roll_list[[h]] <- data.table(
    model = "GLM_QuasiPoisson", cutoff = cutoff,
    train_through = unname(calendar_labels[as.character(cutoff)]),
    test_period = unname(calendar_labels[as.character(cutoff + 1L)]),
    acc = split$test$acc, dev = split$test$dev, cal = split$test$cal,
    actual = split$test$pmts, predicted = ans$prediction
  )
}
qp_rolling_predictions <- rbindlist(qp_roll_list)
stopifnot(nrow(qp_rolling_predictions) == 57L)
qp_metrics <- calculate_metrics(qp_rolling_predictions)
qp_metrics_by_period <- calculate_metrics_by_period(qp_rolling_predictions)

qp_final <- fit_predict_quasipoisson(obs_dat, future_dat)
qp_results <- make_reserve_tables("GLM_QuasiPoisson", qp_final$prediction)
qp_dispersion <- sum(residuals(qp_final$fit, type = "pearson")^2) /
  qp_final$fit$df.residual

print(qp_metrics)
print(qp_results$by_origin)
print(qp_results$totals)
cat("Quasi-Poisson dispersion:", qp_dispersion, "\n")

# ============================================================
# B. TWEEDIE GLM: COMPLETE ANALYSIS

# ============================================================
# TWEEDIE GLM
# Fixed variance power p = 1.38
# obs_dat: observed incremental cells
# future_dat: lower-triangle cells
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(statmod)
})

# ------------------------------------------------------------
# 1. Required-object checks
# ------------------------------------------------------------

required_objects <- c(
  "obs_dat",
  "future_dat",
  "origin_labels",
  "calendar_labels"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

stopifnot(
  nrow(obs_dat) == 105L,
  nrow(future_dat) == 91L,
  all(c("pmts", "acc", "dev", "cal") %in% names(obs_dat)),
  all(c("acc", "dev", "cal") %in% names(future_dat))
)

# ------------------------------------------------------------
# 2. Fixed Tweedie variance power and rolling cutoffs
# ------------------------------------------------------------

tweedie_p_fixed <- 1.38
outer_cutoffs <- 8:13

# ------------------------------------------------------------
# 3. Calendar-based rolling-origin split
# ------------------------------------------------------------

tw_calendar_split <- function(data, cutoff) {
  
  train_data <- data[
    cal <= cutoff
  ]
  
  test_data <- data[
    cal == cutoff + 1L &
      dev > 1L &
      dev <= cutoff
  ]
  
  if (nrow(train_data) == 0L) {
    stop(
      "Empty training sample at cutoff ",
      cutoff
    )
  }
  
  if (nrow(test_data) == 0L) {
    stop(
      "Empty test sample at cutoff ",
      cutoff
    )
  }
  
  list(
    train = train_data,
    test = test_data
  )
}

# ------------------------------------------------------------
# 4. Prepare origin and development factors
# ------------------------------------------------------------

tw_prepare_data <- function(
    train_data,
    test_data
) {
  
  train_copy <- copy(train_data)
  test_copy <- copy(test_data)
  
  train_copy[
    ,
    `:=`(
      acc_f = factor(acc),
      dev_f = factor(dev)
    )
  ]
  
  test_copy[
    ,
    `:=`(
      acc_f = factor(
        acc,
        levels = levels(train_copy$acc_f)
      ),
      dev_f = factor(
        dev,
        levels = levels(train_copy$dev_f)
      )
    )
  ]
  
  if (
    anyNA(test_copy$acc_f) ||
    anyNA(test_copy$dev_f)
  ) {
    stop(
      paste0(
        "The test sample contains an origin or ",
        "development level absent from training."
      )
    )
  }
  
  list(
    train = train_copy,
    test = test_copy
  )
}

# ------------------------------------------------------------
# 5. Fit and predict Tweedie GLM
# ------------------------------------------------------------

tw_fit_predict <- function(
    train_data,
    test_data,
    tweedie_p = 1.38
) {
  
  prepared <- tw_prepare_data(
    train_data = train_data,
    test_data = test_data
  )
  
  fit <- glm(
    pmts ~ acc_f + dev_f,
    family = statmod::tweedie(
      var.power = tweedie_p,
      link.power = 0
    ),
    data = prepared$train,
    control = glm.control(
      epsilon = 1e-8,
      maxit = 500L,
      trace = FALSE
    )
  )
  
  if (!isTRUE(fit$converged)) {
    stop(
      "Tweedie GLM did not converge."
    )
  }
  
  prediction <- predict(
    fit,
    newdata = prepared$test,
    type = "response"
  )
  
  prediction <- pmax(
    as.numeric(prediction),
    0
  )
  
  if (
    length(prediction) != nrow(test_data) ||
    anyNA(prediction) ||
    any(!is.finite(prediction))
  ) {
    stop(
      "Tweedie GLM produced invalid predictions."
    )
  }
  
  list(
    fit = fit,
    prediction = prediction,
    tweedie_p = tweedie_p
  )
}

# ------------------------------------------------------------
# 6. Rolling-origin predictions
# ------------------------------------------------------------

tw_rolling_list <- vector(
  "list",
  length(outer_cutoffs)
)

for (h in seq_along(outer_cutoffs)) {
  
  cutoff <- outer_cutoffs[h]
  
  cat(
    "Tweedie cutoff:",
    cutoff,
    "| train through:",
    unname(
      calendar_labels[
        as.character(cutoff)
      ]
    ),
    "| test period:",
    unname(
      calendar_labels[
        as.character(cutoff + 1L)
      ]
    ),
    "\n"
  )
  
  split_h <- tw_calendar_split(
    data = obs_dat,
    cutoff = cutoff
  )
  
  result_h <- tw_fit_predict(
    train_data = split_h$train,
    test_data = split_h$test,
    tweedie_p = tweedie_p_fixed
  )
  
  tw_rolling_list[[h]] <- data.table(
    model = "GLM_Tweedie",
    cutoff = cutoff,
    
    train_through = unname(
      calendar_labels[
        as.character(cutoff)
      ]
    ),
    
    test_period = unname(
      calendar_labels[
        as.character(cutoff + 1L)
      ]
    ),
    
    acc = split_h$test$acc,
    dev = split_h$test$dev,
    cal = split_h$test$cal,
    actual = split_h$test$pmts,
    predicted = result_h$prediction,
    tweedie_p = tweedie_p_fixed
  )
}

tw_rolling_predictions <- rbindlist(
  tw_rolling_list,
  use.names = TRUE
)

stopifnot(
  nrow(tw_rolling_predictions) == 57L,
  !anyNA(tw_rolling_predictions$actual),
  !anyNA(tw_rolling_predictions$predicted),
  all(is.finite(tw_rolling_predictions$predicted))
)

# ------------------------------------------------------------
# 7. Pooled rolling-origin metrics
# ------------------------------------------------------------

tw_metrics <- tw_rolling_predictions[
  ,
  {
    SSE <- sum(
      (actual - predicted)^2
    )
    
    SST <- sum(
      (actual - mean(actual))^2
    )
    
    .(
      N_test = .N,
      
      RMSE = sqrt(
        mean(
          (actual - predicted)^2
        )
      ),
      
      MAE = mean(
        abs(actual - predicted)
      ),
      
      R_squared_predictive =
        1 - SSE / SST,
      
      Actual_total =
        sum(actual),
      
      Predicted_total =
        sum(predicted),
      
      Aggregate_relative_error =
        (
          sum(predicted) -
            sum(actual)
        ) /
        sum(actual)
    )
  }
]

cat(
  "\n--- TWEEDIE POOLED ROLLING-ORIGIN METRICS ---\n"
)

print(tw_metrics)

# ------------------------------------------------------------
# 8. Metrics by calendar test period
# ------------------------------------------------------------

tw_metrics_by_period <-
  tw_rolling_predictions[
    ,
    .(
      N_test = .N,
      
      RMSE = sqrt(
        mean(
          (actual - predicted)^2
        )
      ),
      
      MAE = mean(
        abs(actual - predicted)
      ),
      
      Actual_total =
        sum(actual),
      
      Predicted_total =
        sum(predicted),
      
      Aggregate_relative_error =
        (
          sum(predicted) -
            sum(actual)
        ) /
        sum(actual)
    ),
    by = .(
      cutoff,
      train_through,
      test_period
    )
  ]

cat(
  "\n--- TWEEDIE METRICS BY CALENDAR PERIOD ---\n"
)

print(tw_metrics_by_period)

# ------------------------------------------------------------
# 9. Final fit on all 105 observed cells
# ------------------------------------------------------------

tw_final <- tw_fit_predict(
  train_data = obs_dat,
  test_data = future_dat,
  tweedie_p = tweedie_p_fixed
)

tw_future_predictions <- cbind(
  copy(future_dat),
  data.table(
    predicted_incremental =
      tw_final$prediction
  )
)

tw_future_predictions[
  ,
  `:=`(
    model = "GLM_Tweedie",
    origin_period =
      origin_labels[acc]
  )
]

setorder(
  tw_future_predictions,
  acc,
  dev
)

# ------------------------------------------------------------
# 10. Latest observed cumulative paid claims
# ------------------------------------------------------------

tw_latest_observed <- obs_dat[
  ,
  .(
    observed_cumulative_paid =
      sum(pmts)
  ),
  by = acc
]

tw_latest_observed[
  ,
  origin_period :=
    origin_labels[acc]
]

# ------------------------------------------------------------
# 11. Reserve by origin period
# ------------------------------------------------------------

tw_reserve_by_origin <-
  tw_future_predictions[
    ,
    .(
      reserve =
        sum(predicted_incremental),
      
      future_cells = .N
    ),
    by = .(
      model,
      acc,
      origin_period
    )
  ]

tw_all_origins <- data.table(
  model = "GLM_Tweedie",
  acc = seq_along(origin_labels),
  origin_period = origin_labels
)

tw_reserve_by_origin <- merge(
  tw_all_origins,
  tw_reserve_by_origin,
  by = c(
    "model",
    "acc",
    "origin_period"
  ),
  all.x = TRUE
)

tw_reserve_by_origin[
  is.na(reserve),
  `:=`(
    reserve = 0,
    future_cells = 0L
  )
]

tw_reserve_by_origin <- merge(
  tw_reserve_by_origin,
  tw_latest_observed,
  by = c(
    "acc",
    "origin_period"
  ),
  all.x = TRUE
)

tw_reserve_by_origin[
  ,
  ultimate_claims :=
    observed_cumulative_paid +
    reserve
]

setorder(
  tw_reserve_by_origin,
  acc
)

cat(
  "\n--- TWEEDIE RESERVE AND ULTIMATE CLAIMS BY ORIGIN ---\n"
)

print(tw_reserve_by_origin)

# ------------------------------------------------------------
# 12. Total reserve and ultimate claims
# ------------------------------------------------------------

tw_reserve_totals <-
  tw_reserve_by_origin[
    ,
    .(
      observed_cumulative_paid =
        sum(observed_cumulative_paid),
      
      reserve =
        sum(reserve),
      
      ultimate_claims =
        sum(ultimate_claims)
    ),
    by = model
  ]

cat(
  "\n--- TWEEDIE TOTAL RESERVE AND ULTIMATE CLAIMS ---\n"
)

print(tw_reserve_totals)



# ============================================================
# C.LASSO-penalized Poisson regression for incremental paid claims; alpha = 1 gives the L1 penalty
#===================================================================

LinearSpline <- function(x, start, stop) {
  pmin(stop - start, pmax(0, x - start))
}

GetScaling <- function(x) {
  value <- sqrt(sum((x - mean(x))^2) / length(x))
  if (!is.finite(value) || value <= 0) 1 else value
}

GetRamps <- function(x, name, np, scaling) {
  out <- matrix(NA_real_, nrow = length(x), ncol = np - 1L)
  for (i in seq_len(np - 1L)) out[, i] <- LinearSpline(x, i, 999) / scaling
  colnames(out) <- paste0("L_", seq_len(np - 1L), "_999_", name)
  out
}

GetInts <- function(x1, x2, name1, name2, np, scale1, scale2) {
  out <- matrix(NA_real_, nrow = length(x1), ncol = (np - 1L)^2)
  names_out <- character((np - 1L)^2)
  k <- 0L
  for (i in 2:np) {
    a <- LinearSpline(x1, i - 1L, i) / scale1
    for (j in 2:np) {
      k <- k + 1L
      b <- LinearSpline(x2, j - 1L, j) / scale2
      out[, k] <- a * b
      names_out[k] <- paste0(name1, "_ge_", i, "_x_", name2, "_ge_", j)
    }
  }
  colnames(out) <- names_out
  out
}

make_lasso_xy <- function(train, test) {
  s_acc <- GetScaling(train$acc)
  s_dev <- GetScaling(train$dev)
  s_cal <- GetScaling(train$cal)

  build_x <- function(dat) {
    x <- cbind(
      GetRamps(dat$acc, "acc", N, s_acc),
      GetRamps(dat$dev, "dev", N, s_dev),
      GetRamps(dat$cal, "cal", N, s_cal),
      GetInts(dat$acc, dat$dev, "acc", "dev", N, s_acc, s_dev),
      GetInts(dat$dev, dat$cal, "dev", "cal", N, s_dev, s_cal),
      GetInts(dat$acc, dat$cal, "acc", "cal", N, s_acc, s_cal)
    )
    storage.mode(x) <- "double"
    x
  }

  x_train <- build_x(train)
  x_test <- build_x(test)
  train_sd <- apply(x_train, 2L, sd)
  keep <- is.finite(train_sd) & train_sd > 0
  if (!any(keep)) stop("No non-constant LASSO features remain.")
  list(
    x_train = x_train[, keep, drop = FALSE],
    x_test = x_test[, keep, drop = FALSE]
  )
}

lasso_grid <- data.table(
  lambda_ratio = c(0.001, 0.003, 0.010, 0.030, 0.100)
)
lasso_grid[, candidate_id := .I]

fit_predict_lasso <- function(train, test, lambda_ratio) {
  xy <- make_lasso_xy(train, test)
  initial_path <- glmnet(
    x = xy$x_train, y = train$pmts,
    family = "poisson", alpha = 1, standardize = FALSE,
    intercept = TRUE, nlambda = 5L, lambda.min.ratio = 0.10,
    maxit = 200000
  )
  lambda_max <- max(initial_path$lambda)
  lambda_value <- lambda_max * lambda_ratio
  lambda_sequence <- sort(unique(c(lambda_max, lambda_value)), decreasing = TRUE)

  fit <- glmnet(
    x = xy$x_train, y = train$pmts,
    family = "poisson", alpha = 1,
    lambda = lambda_sequence, standardize = FALSE,
    intercept = TRUE, thresh = 1e-8, maxit = 200000
  )
  pred <- pmax(as.numeric(predict(
    fit, newx = xy$x_test, s = lambda_value, type = "response"
  )), 0)
  if (length(pred) != nrow(test) || anyNA(pred) || any(!is.finite(pred))) {
    stop("LASSO produced invalid predictions.")
  }
  list(
    prediction = pred, lambda = lambda_value,
    lambda_max = lambda_max, lambda_ratio = lambda_ratio
  )
}

# Inner tuning is performed using only the corresponding outer training data.
tune_lasso_inner <- function(outer_cutoff) {
  inner_cutoffs <- inner_start_cutoff:(outer_cutoff - 1L)
  scores <- vector("list", nrow(lasso_grid))
  for (g in seq_len(nrow(lasso_grid))) {
    fold_predictions <- vector("list", length(inner_cutoffs))
    for (h in seq_along(inner_cutoffs)) {
      split <- calendar_split(obs_dat[cal <= outer_cutoff], inner_cutoffs[h])
      ans <- fit_predict_lasso(
        split$train, split$test, lasso_grid$lambda_ratio[g]
      )
      fold_predictions[[h]] <- data.table(
        actual = split$test$pmts, predicted = ans$prediction
      )
    }
    pp <- rbindlist(fold_predictions)
    scores[[g]] <- data.table(
      candidate_id = lasso_grid$candidate_id[g],
      inner_N = nrow(pp),
      inner_RMSE = sqrt(mean((pp$actual - pp$predicted)^2)),
      inner_MAE = mean(abs(pp$actual - pp$predicted))
    )
  }
  ans <- merge(lasso_grid, rbindlist(scores), by = "candidate_id")
  setorder(ans, inner_RMSE, inner_MAE, candidate_id)
  ans
}

lasso_roll_list <- vector("list", length(outer_cutoffs))
lasso_best_outer_list <- vector("list", length(outer_cutoffs))
lasso_all_inner_list <- vector("list", length(outer_cutoffs))

for (h in seq_along(outer_cutoffs)) {
  cutoff <- outer_cutoffs[h]
  split <- calendar_split(obs_dat, cutoff)
  tuning <- tune_lasso_inner(cutoff)
  best <- tuning[1L]
  ans <- fit_predict_lasso(split$train, split$test, best$lambda_ratio)

  lasso_roll_list[[h]] <- data.table(
    model = "LASSO_L1", cutoff = cutoff,
    train_through = unname(calendar_labels[as.character(cutoff)]),
    test_period = unname(calendar_labels[as.character(cutoff + 1L)]),
    acc = split$test$acc, dev = split$test$dev, cal = split$test$cal,
    actual = split$test$pmts, predicted = ans$prediction
  )
  lasso_best_outer_list[[h]] <- cbind(
    data.table(outer_cutoff = cutoff), best,
    data.table(selected_lambda = ans$lambda)
  )
  lasso_all_inner_list[[h]] <- cbind(
    data.table(outer_cutoff = cutoff), tuning
  )
}

lasso_rolling_predictions <- rbindlist(lasso_roll_list)
lasso_best_by_outer_fold <- rbindlist(lasso_best_outer_list)
lasso_all_inner_scores <- rbindlist(lasso_all_inner_list)
stopifnot(nrow(lasso_rolling_predictions) == 57L)
lasso_metrics <- calculate_metrics(lasso_rolling_predictions)
lasso_metrics_by_period <- calculate_metrics_by_period(lasso_rolling_predictions)

# Select the final lambda ratio from all six historical validation diagonals.
lasso_final_score_list <- vector("list", nrow(lasso_grid))
for (g in seq_len(nrow(lasso_grid))) {
  fold_predictions <- vector("list", length(outer_cutoffs))
  for (h in seq_along(outer_cutoffs)) {
    split <- calendar_split(obs_dat, outer_cutoffs[h])
    ans <- fit_predict_lasso(
      split$train, split$test, lasso_grid$lambda_ratio[g]
    )
    fold_predictions[[h]] <- data.table(
      actual = split$test$pmts, predicted = ans$prediction
    )
  }
  pp <- rbindlist(fold_predictions)
  lasso_final_score_list[[g]] <- data.table(
    candidate_id = lasso_grid$candidate_id[g],
    tuning_N = nrow(pp),
    tuning_RMSE = sqrt(mean((pp$actual - pp$predicted)^2)),
    tuning_MAE = mean(abs(pp$actual - pp$predicted))
  )
}
lasso_final_scores <- merge(
  lasso_grid, rbindlist(lasso_final_score_list), by = "candidate_id"
)
setorder(lasso_final_scores, tuning_RMSE, tuning_MAE, candidate_id)
lasso_final_best <- lasso_final_scores[1L]

lasso_final_fit <- fit_predict_lasso(
  obs_dat, future_dat, lasso_final_best$lambda_ratio
)
lasso_results <- make_reserve_tables("LASSO_L1", lasso_final_fit$prediction)
lasso_final_parameters <- cbind(
  lasso_final_best,
  data.table(
    alpha = 1,
    final_lambda_max = lasso_final_fit$lambda_max,
    final_lambda = lasso_final_fit$lambda
  )
)

print(lasso_metrics)
print(lasso_final_parameters)
print(lasso_results$by_origin)
print(lasso_results$totals)

############    END   #############

