# ============================================================
# MACHINE-LEARNING CLAIM RESERVING
# Random Forest -> Decision Tree -> XGBoost
#
# Target: incremental paid claims
# Predictors: origin (acc), development (dev), calendar (cal)
# Validation: nested rolling-origin by calendar period
# Final reserve: model fitted on all 105 observed cells
# ============================================================


suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(rpart)
  library(xgboost)
})

# Reproducibility. The models are run with one processing thread.
set.seed(123)

# ============================================================
# 1. CUMULATIVE PAID-CLAIMS TRIANGLE
# ============================================================

triangle_data <- matrix(
  c(
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
  ),
  nrow = 14L,
  byrow = TRUE
)

N <- nrow(triangle_data)
stopifnot(N == 14L, ncol(triangle_data) == N)

origin_labels <- c(
  "2018_S1", "2018_S2", "2019_S1", "2019_S2",
  "2020_S1", "2020_S2", "2021_S1", "2021_S2",
  "2022_S1", "2022_S2", "2023_S1", "2023_S2",
  "2024_S1", "2024_S2"
)

calendar_labels <- setNames(
  c(
    "2018_S1", "2018_S2", "2019_S1", "2019_S2",
    "2020_S1", "2020_S2", "2021_S1", "2021_S2",
    "2022_S1", "2022_S2", "2023_S1", "2023_S2",
    "2024_S1", "2024_S2"
  ),
  as.character(1:14)
)

# ============================================================
# 2. CUMULATIVE TO INCREMENTAL; WIDE TO LONG
# ============================================================

incremental_triangle <- triangle_data
for (i in seq_len(N)) {
  for (j in seq_len(N)) {
    if (is.na(triangle_data[i, j])) {
      incremental_triangle[i, j] <- NA_real_
    } else if (j == 1L) {
      incremental_triangle[i, j] <- triangle_data[i, j]
    } else {
      incremental_triangle[i, j] <-
        triangle_data[i, j] - triangle_data[i, j - 1L]
    }
  }
}

long_dat <- CJ(acc = seq_len(N), dev = seq_len(N))
long_dat[, pmts := as.numeric(incremental_triangle[cbind(acc, dev)])]
long_dat[, cal := acc + dev - 1L]
long_dat[, observed := !is.na(pmts)]
setorder(long_dat, acc, dev)

obs_dat <- long_dat[observed == TRUE, .(acc, dev, cal, pmts)]
future_dat <- long_dat[observed == FALSE, .(acc, dev, cal)]

stopifnot(
  nrow(obs_dat) == 105L,
  nrow(future_dat) == 91L,
  all(obs_dat$pmts >= 0),
  !anyNA(obs_dat)
)

cat("Observed cells:", nrow(obs_dat), "\n")
cat("Future cells:", nrow(future_dat), "\n")

# ============================================================
# 3. COMMON ROLLING-ORIGIN DEFINITIONS AND HELPERS
# ============================================================

# These outer test periods are identical for RF, DT and XGBoost.
outer_cutoffs <- 8:13

calendar_split <- function(data, cutoff) {
  train <- data[cal <= cutoff]

  # dev > 1 excludes the first development cell, whose origin level
  # is not represented in the historical training triangle.
  # dev <= cutoff creates the common 57-cell comparison sample.
  test <- data[
    cal == cutoff + 1L &
      dev > 1L &
      dev <= cutoff
  ]

  stopifnot(nrow(train) > 0L, nrow(test) > 0L)
  list(train = train, test = test)
}

# Inner historical cutoffs used inside a given outer training sample.
inner_cutoffs_for <- function(outer_cutoff) {
  cutoffs <- 5:(outer_cutoff - 1L)
  cutoffs[cutoffs >= 5L]
}

predictive_metrics <- function(predictions) {
  stopifnot(
    nrow(predictions) > 0L,
    !anyNA(predictions$actual),
    !anyNA(predictions$predicted),
    all(is.finite(predictions$predicted))
  )

  sse <- sum((predictions$actual - predictions$predicted)^2)
  sst <- sum((predictions$actual - mean(predictions$actual))^2)

  data.table(
    N_test = nrow(predictions),
    RMSE = sqrt(mean((predictions$actual - predictions$predicted)^2)),
    MAE = mean(abs(predictions$actual - predictions$predicted)),
    R_squared_predictive = 1 - sse / sst,
    Actual_total = sum(predictions$actual),
    Predicted_total = sum(predictions$predicted),
    Aggregate_relative_error =
      (sum(predictions$predicted) - sum(predictions$actual)) /
      sum(predictions$actual)
  )
}

metrics_by_period <- function(predictions) {
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

make_reserve_tables <- function(model_name, future_predictions) {
  stopifnot(
    length(future_predictions) == nrow(future_dat),
    !anyNA(future_predictions),
    all(is.finite(future_predictions)),
    all(future_predictions >= 0)
  )

  cell_predictions <- cbind(
    copy(future_dat),
    data.table(predicted_incremental = future_predictions)
  )
  cell_predictions[, `:=`(
    model = model_name,
    origin_period = origin_labels[acc]
  )]
  setorder(cell_predictions, acc, dev)

  observed_paid <- obs_dat[, .(
    observed_cumulative_paid = sum(pmts)
  ), by = acc]

  by_origin <- cell_predictions[, .(
    reserve = sum(predicted_incremental),
    future_cells = .N
  ), by = .(model, acc, origin_period)]

  all_origins <- data.table(
    model = model_name,
    acc = seq_len(N),
    origin_period = origin_labels
  )

  by_origin <- merge(
    all_origins, by_origin,
    by = c("model", "acc", "origin_period"),
    all.x = TRUE
  )
  by_origin[is.na(reserve), `:=`(reserve = 0, future_cells = 0L)]
  by_origin <- merge(by_origin, observed_paid, by = "acc", all.x = TRUE)
  by_origin[, ultimate_claims := observed_cumulative_paid + reserve]
  setorder(by_origin, acc)

  totals <- by_origin[, .(
    observed_cumulative_paid = sum(observed_cumulative_paid),
    reserve = sum(reserve),
    ultimate_claims = sum(ultimate_claims)
  ), by = model]

  list(
    cell_predictions = cell_predictions,
    by_origin = by_origin,
    totals = totals
  )
}

# Check the common outer test sample before fitting any model.
common_test_sample <- rbindlist(lapply(outer_cutoffs, function(cutoff) {
  split <- calendar_split(obs_dat, cutoff)
  data.table(
    cutoff = cutoff,
    train_through = unname(calendar_labels[as.character(cutoff)]),
    test_period = unname(calendar_labels[as.character(cutoff + 1L)]),
    N_test = nrow(split$test),
    Actual_total = sum(split$test$pmts)
  )
}))

stopifnot(
  sum(common_test_sample$N_test) == 57L,
  sum(common_test_sample$Actual_total) == 193548009
)

cat("\n--- COMMON OUTER ROLLING-ORIGIN TEST SAMPLE ---\n")
print(common_test_sample)

# ============================================================
# 4. RANDOM FOREST: TUNING, METRICS AND FINAL RESERVE
# ============================================================

rf_grid <- CJ(
  mtry = 1:3,
  min.node.size = c(3L, 5L, 8L),
  sample.fraction = c(0.7, 0.9, 1.0)
)
rf_grid[, `:=`(candidate_id = .I, num.trees = 500L)]
setcolorder(rf_grid, c(
  "candidate_id", "mtry", "min.node.size",
  "sample.fraction", "num.trees"
))

rf_fit_predict <- function(train_data, test_data, pars) {
  fit <- ranger(
    pmts ~ acc + dev + cal,
    data = as.data.frame(train_data),
    num.trees = pars$num.trees,
    mtry = pars$mtry,
    min.node.size = pars$min.node.size,
    sample.fraction = pars$sample.fraction,
    replace = TRUE,
    seed = 123,
    num.threads = 1,
    write.forest = TRUE
  )
  pmax(as.numeric(predict(fit, data = as.data.frame(test_data))$predictions), 0)
}

rf_tune <- function(training_data, inner_cutoffs) {
  candidate_scores <- vector("list", nrow(rf_grid))

  for (g in seq_len(nrow(rf_grid))) {
    fold_predictions <- vector("list", length(inner_cutoffs))

    for (h in seq_along(inner_cutoffs)) {
      split <- calendar_split(training_data, inner_cutoffs[h])
      pred <- rf_fit_predict(split$train, split$test, rf_grid[g])
      fold_predictions[[h]] <- data.table(
        actual = split$test$pmts,
        predicted = pred
      )
    }

    pooled <- rbindlist(fold_predictions)
    candidate_scores[[g]] <- data.table(
      candidate_id = rf_grid$candidate_id[g],
      tuning_N = nrow(pooled),
      tuning_RMSE = sqrt(mean((pooled$actual - pooled$predicted)^2)),
      tuning_MAE = mean(abs(pooled$actual - pooled$predicted))
    )
  }

  scores <- merge(rf_grid, rbindlist(candidate_scores), by = "candidate_id")
  setorder(scores, tuning_RMSE, tuning_MAE, candidate_id)
  list(best = scores[1L], scores = scores)
}

# outer performance: tuning uses only the history available
# inside each outer training sample.
rf_outer_list <- vector("list", length(outer_cutoffs))
rf_outer_best_list <- vector("list", length(outer_cutoffs))

for (h in seq_along(outer_cutoffs)) {
  cutoff <- outer_cutoffs[h]
  outer_split <- calendar_split(obs_dat, cutoff)
  tuned <- rf_tune(outer_split$train, inner_cutoffs_for(cutoff))
  best <- tuned$best

  pred <- rf_fit_predict(outer_split$train, outer_split$test, best)

  rf_outer_list[[h]] <- data.table(
    model = "RF",
    cutoff = cutoff,
    train_through = unname(calendar_labels[as.character(cutoff)]),
    test_period = unname(calendar_labels[as.character(cutoff + 1L)]),
    acc = outer_split$test$acc,
    dev = outer_split$test$dev,
    cal = outer_split$test$cal,
    actual = outer_split$test$pmts,
    predicted = pred
  )

  rf_outer_best_list[[h]] <- cbind(
    data.table(outer_cutoff = cutoff), best
  )
}

rf_rolling_predictions <- rbindlist(rf_outer_list)
rf_outer_best_parameters <- rbindlist(rf_outer_best_list)
rf_metrics <- predictive_metrics(rf_rolling_predictions)
rf_metrics_period <- metrics_by_period(rf_rolling_predictions)

stopifnot(nrow(rf_rolling_predictions) == 57L)

# Final parameter selection uses all six historical rolling folds.
rf_final_tuning <- rf_tune(obs_dat, outer_cutoffs)
rf_final_parameters <- rf_final_tuning$best
rf_final_prediction <- rf_fit_predict(
  obs_dat, future_dat, rf_final_parameters
)
rf_results <- make_reserve_tables("RF", rf_final_prediction)

cat("\n--- RF FINAL OPTIMAL PARAMETERS ---\n")
print(rf_final_parameters)
cat("\n--- RF HONEST ROLLING-ORIGIN METRICS ---\n")
print(rf_metrics)
cat("\n--- RF METRICS BY CALENDAR PERIOD ---\n")
print(rf_metrics_period)
cat("\n--- RF RESERVE AND ULTIMATE CLAIMS BY ORIGIN ---\n")
print(rf_results$by_origin)
cat("\n--- RF TOTAL RESERVE AND ULTIMATE CLAIMS ---\n")
print(rf_results$totals)


# ============================================================
# 5. DECISION TREE: TUNING, METRICS AND FINAL RESERVE
# ============================================================

dt_grid <- CJ(
  cp = c(0.001, 0.005, 0.010),
  minsplit = c(4L, 6L)
)
dt_grid[, candidate_id := .I]
setcolorder(dt_grid, c("candidate_id", "cp", "minsplit"))

dt_fit_predict <- function(train_data, test_data, pars) {
  fit <- rpart(
    pmts ~ acc + dev + cal,
    data = as.data.frame(train_data),
    method = "anova",
    control = rpart.control(
      cp = pars$cp,
      minsplit = pars$minsplit,
      xval = 0L
    )
  )
  pmax(as.numeric(predict(fit, newdata = as.data.frame(test_data))), 0)
}

dt_tune <- function(training_data, inner_cutoffs) {
  candidate_scores <- vector("list", nrow(dt_grid))

  for (g in seq_len(nrow(dt_grid))) {
    fold_predictions <- vector("list", length(inner_cutoffs))

    for (h in seq_along(inner_cutoffs)) {
      split <- calendar_split(training_data, inner_cutoffs[h])
      pred <- dt_fit_predict(split$train, split$test, dt_grid[g])
      fold_predictions[[h]] <- data.table(
        actual = split$test$pmts,
        predicted = pred
      )
    }

    pooled <- rbindlist(fold_predictions)
    candidate_scores[[g]] <- data.table(
      candidate_id = dt_grid$candidate_id[g],
      tuning_N = nrow(pooled),
      tuning_RMSE = sqrt(mean((pooled$actual - pooled$predicted)^2)),
      tuning_MAE = mean(abs(pooled$actual - pooled$predicted))
    )
  }

  scores <- merge(dt_grid, rbindlist(candidate_scores), by = "candidate_id")
  setorder(scores, tuning_RMSE, tuning_MAE, candidate_id)
  list(best = scores[1L], scores = scores)
}

dt_outer_list <- vector("list", length(outer_cutoffs))
dt_outer_best_list <- vector("list", length(outer_cutoffs))

for (h in seq_along(outer_cutoffs)) {
  cutoff <- outer_cutoffs[h]
  outer_split <- calendar_split(obs_dat, cutoff)
  tuned <- dt_tune(outer_split$train, inner_cutoffs_for(cutoff))
  best <- tuned$best
  pred <- dt_fit_predict(outer_split$train, outer_split$test, best)

  dt_outer_list[[h]] <- data.table(
    model = "DT",
    cutoff = cutoff,
    train_through = unname(calendar_labels[as.character(cutoff)]),
    test_period = unname(calendar_labels[as.character(cutoff + 1L)]),
    acc = outer_split$test$acc,
    dev = outer_split$test$dev,
    cal = outer_split$test$cal,
    actual = outer_split$test$pmts,
    predicted = pred
  )

  dt_outer_best_list[[h]] <- cbind(
    data.table(outer_cutoff = cutoff), best
  )
}

dt_rolling_predictions <- rbindlist(dt_outer_list)
dt_outer_best_parameters <- rbindlist(dt_outer_best_list)
dt_metrics <- predictive_metrics(dt_rolling_predictions)
dt_metrics_period <- metrics_by_period(dt_rolling_predictions)

stopifnot(nrow(dt_rolling_predictions) == 57L)

dt_final_tuning <- dt_tune(obs_dat, outer_cutoffs)
dt_final_parameters <- dt_final_tuning$best
dt_final_prediction <- dt_fit_predict(
  obs_dat, future_dat, dt_final_parameters
)
dt_results <- make_reserve_tables("DT", dt_final_prediction)

cat("\n--- DT FINAL OPTIMAL PARAMETERS ---\n")
print(dt_final_parameters)
cat("\n--- DT HONEST ROLLING-ORIGIN METRICS ---\n")
print(dt_metrics)
cat("\n--- DT METRICS BY CALENDAR PERIOD ---\n")
print(dt_metrics_period)
cat("\n--- DT RESERVE AND ULTIMATE CLAIMS BY ORIGIN ---\n")
print(dt_results$by_origin)
cat("\n--- DT TOTAL RESERVE AND ULTIMATE CLAIMS ---\n")
print(dt_results$totals)

########################

# ============================================================
# XGBOOST
# Tuning -> nested rolling-origin evaluation -> final reserve
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
})

# ------------------------------------------------------------
# 1. XGBoost tuning grid
# ------------------------------------------------------------

xgb_grid <- CJ(
  nrounds = c(300L, 600L),
  eta = c(0.05, 0.10),
  max_depth = c(2L, 3L),
  min_child_weight = c(1, 5),
  subsample = c(0.7, 1.0)
)

xgb_grid[
  ,
  `:=`(
    candidate_id = .I,
    colsample_bytree = 1.0,
    gamma = 0.0
  )
]

setcolorder(
  xgb_grid,
  c(
    "candidate_id",
    "nrounds",
    "eta",
    "max_depth",
    "min_child_weight",
    "subsample",
    "colsample_bytree",
    "gamma"
  )
)

cat("\n--- XGBOOST TUNING GRID ---\n")
print(xgb_grid)

# ------------------------------------------------------------
# 2. Convert predictors to an XGBoost matrix
# ------------------------------------------------------------

xgb_matrix <- function(data) {
  
  matrix_x <- as.matrix(
    data[, .(acc, dev, cal)]
  )
  
  storage.mode(matrix_x) <- "double"
  
  matrix_x
}

# ------------------------------------------------------------
# 3. Fit XGBoost and generate predictions
# ------------------------------------------------------------

xgb_fit_predict <- function(
    train_data,
    test_data,
    parameters
) {
  
  stopifnot(
    nrow(train_data) > 0L,
    nrow(test_data) > 0L,
    !anyNA(train_data[, .(acc, dev, cal, pmts)]),
    !anyNA(test_data[, .(acc, dev, cal)])
  )
  
  dtrain <- xgb.DMatrix(
    data = xgb_matrix(train_data),
    label = train_data$pmts
  )
  
  dtest <- xgb.DMatrix(
    data = xgb_matrix(test_data)
  )
  
  # XGBoost reproducibility is controlled from R.

  set.seed(123)
  
  fitted_model <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      
      eta = as.numeric(
        parameters$eta
      ),
      
      max_depth = as.integer(
        parameters$max_depth
      ),
      
      min_child_weight = as.numeric(
        parameters$min_child_weight
      ),
      
      subsample = as.numeric(
        parameters$subsample
      ),
      
      colsample_bytree = as.numeric(
        parameters$colsample_bytree
      ),

      gamma = as.numeric(
        parameters$gamma
      ),
      
      nthread = 1
    ),
    
    data = dtrain,
    
    nrounds = as.integer(
      parameters$nrounds
    ),
    
    verbose = 0
  )
  
  predictions <- as.numeric(
    predict(
      fitted_model,
      dtest
    )
  )
  
  # Incremental paid claims cannot be negative.
  predictions <- pmax(
    predictions,
    0
  )
  
  if (
    length(predictions) != nrow(test_data) ||
    anyNA(predictions) ||
    any(!is.finite(predictions))
  ) {
    stop(
      "XGBoost produced invalid predictions."
    )
  }
  
  predictions
}

# ------------------------------------------------------------
# 4. Tune XGBoost using historical rolling-origin folds
# ------------------------------------------------------------

xgb_tune <- function(
    training_data,
    inner_cutoffs
) {
  
  candidate_results <- vector(
    mode = "list",
    length = nrow(xgb_grid)
  )
  
  for (g in seq_len(nrow(xgb_grid))) {
    
    fold_predictions <- vector(
      mode = "list",
      length = length(inner_cutoffs)
    )
    
    for (h in seq_along(inner_cutoffs)) {
      
      cutoff <- inner_cutoffs[h]
      
      split_h <- calendar_split(
        data = training_data,
        cutoff = cutoff
      )
      
      predicted_h <- xgb_fit_predict(
        train_data = split_h$train,
        test_data = split_h$test,
        parameters = xgb_grid[g]
      )
      
      fold_predictions[[h]] <- data.table(
        cutoff = cutoff,
        actual = split_h$test$pmts,
        predicted = predicted_h
      )
    }
    
    pooled_predictions <- rbindlist(
      fold_predictions
    )
    
    candidate_results[[g]] <- data.table(
      candidate_id =
        xgb_grid$candidate_id[g],
      
      tuning_N =
        nrow(pooled_predictions),
      
      tuning_RMSE =
        sqrt(
          mean(
            (
              pooled_predictions$actual -
                pooled_predictions$predicted
            )^2
          )
        ),
      
      tuning_MAE =
        mean(
          abs(
            pooled_predictions$actual -
              pooled_predictions$predicted
          )
        )
    )
  }
  
  tuning_results <- merge(
    xgb_grid,
    rbindlist(candidate_results),
    by = "candidate_id"
  )
  
  setorder(
    tuning_results,
    tuning_RMSE,
    tuning_MAE,
    candidate_id
  )
  
  list(
    best = tuning_results[1L],
    scores = tuning_results
  )
}

# ------------------------------------------------------------
# 5.nested rolling-origin evaluation
# ------------------------------------------------------------

xgb_outer_predictions_list <- vector(
  mode = "list",
  length = length(outer_cutoffs)
)

xgb_outer_parameters_list <- vector(
  mode = "list",
  length = length(outer_cutoffs)
)

for (h in seq_along(outer_cutoffs)) {
  
  outer_cutoff <- outer_cutoffs[h]
  
  cat(
    "\nXGBoost outer cutoff:",
    outer_cutoff,
    "\n"
  )
  
  outer_split <- calendar_split(
    data = obs_dat,
    cutoff = outer_cutoff
  )
  
  # Tuning is conducted only within the outer training sample.
  inner_cutoffs <- inner_cutoffs_for(
    outer_cutoff
  )
  
  xgb_inner_tuning <- xgb_tune(
    training_data = outer_split$train,
    inner_cutoffs = inner_cutoffs
  )
  
  best_parameters <- xgb_inner_tuning$best
  
  outer_prediction <- xgb_fit_predict(
    train_data = outer_split$train,
    test_data = outer_split$test,
    parameters = best_parameters
  )
  
  xgb_outer_predictions_list[[h]] <- data.table(
    model = "XGBoost",
    
    cutoff = outer_cutoff,
    
    train_through =
      unname(
        calendar_labels[
          as.character(outer_cutoff)
        ]
      ),
    
    test_period =
      unname(
        calendar_labels[
          as.character(outer_cutoff + 1L)
        ]
      ),
    
    acc = outer_split$test$acc,
    dev = outer_split$test$dev,
    cal = outer_split$test$cal,
    actual = outer_split$test$pmts,
    predicted = outer_prediction
  )
  
  xgb_outer_parameters_list[[h]] <- cbind(
    data.table(
      outer_cutoff = outer_cutoff
    ),
    best_parameters
  )
}

# ------------------------------------------------------------
# 6. Combine the outer test predictions
# ------------------------------------------------------------

xgb_rolling_predictions <- rbindlist(
  xgb_outer_predictions_list
)

xgb_outer_best_parameters <- rbindlist(
  xgb_outer_parameters_list
)

stopifnot(
  nrow(xgb_rolling_predictions) == 57L,
  
  !anyNA(
    xgb_rolling_predictions$predicted
  ),
  
  all(
    is.finite(
      xgb_rolling_predictions$predicted
    )
  ),
  
  sum(
    xgb_rolling_predictions$actual
  ) == 193548009
)

# ------------------------------------------------------------
# 7. Calculate predictive metrics
# ------------------------------------------------------------

xgb_metrics <- predictive_metrics(
  xgb_rolling_predictions
)

xgb_metrics_period <- metrics_by_period(
  xgb_rolling_predictions
)

cat(
  "\n--- XGBOOST HONEST ROLLING-ORIGIN METRICS ---\n"
)

print(
  xgb_metrics
)

cat(
  "\n--- XGBOOST METRICS BY CALENDAR PERIOD ---\n"
)

print(
  xgb_metrics_period
)

cat(
  "\n--- XGBOOST PARAMETERS SELECTED IN EACH OUTER FOLD ---\n"
)

print(
  xgb_outer_best_parameters
)

# ------------------------------------------------------------
# 8. Final tuning using all historical rolling-origin folds
# ------------------------------------------------------------

xgb_final_tuning <- xgb_tune(
  training_data = obs_dat,
  inner_cutoffs = outer_cutoffs
)

xgb_final_parameters <- xgb_final_tuning$best

cat(
  "\n--- XGBOOST FINAL OPTIMAL PARAMETERS ---\n"
)

print(
  xgb_final_parameters
)

# ------------------------------------------------------------
# 9. Final fit on all 105 observed cells
# ------------------------------------------------------------

xgb_final_prediction <- xgb_fit_predict(
  train_data = obs_dat,
  test_data = future_dat,
  parameters = xgb_final_parameters
)

stopifnot(
  length(xgb_final_prediction) ==
    nrow(future_dat),
  
  !anyNA(xgb_final_prediction),
  
  all(
    is.finite(
      xgb_final_prediction
    )
  ),
  
  all(xgb_final_prediction >= 0)
)

# ------------------------------------------------------------
# 10. Reserve and ultimate claims
# ------------------------------------------------------------

xgb_results <- make_reserve_tables(
  model_name = "XGBoost",
  future_predictions =
    xgb_final_prediction
)

cat(
  "\n--- XGBOOST RESERVE AND ULTIMATE CLAIMS BY ORIGIN ---\n"
)

print(
  xgb_results$by_origin
)

cat(
  "\n--- XGBOOST TOTAL RESERVE AND ULTIMATE CLAIMS ---\n"
)

print(
  xgb_results$totals
)

# ------------------------------------------------------------
# 11. Final consistency checks
# ------------------------------------------------------------

stopifnot(
  nrow(xgb_results$by_origin) == 14L,
  
  nrow(xgb_results$totals) == 1L,
  
  abs(
    sum(xgb_results$by_origin$reserve) -
      xgb_results$totals$reserve
  ) < 0.01,
  
  abs(
    xgb_results$totals$observed_cumulative_paid +
      xgb_results$totals$reserve -
      xgb_results$totals$ultimate_claims
  ) < 0.01
)


##################  END ###################################
