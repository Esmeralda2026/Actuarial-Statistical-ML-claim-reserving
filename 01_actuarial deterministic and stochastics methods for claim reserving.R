
#-------------- Actuarial Methods for claim reserving--------------#

#---------------Deterministic methods------------#
######################################################
#---------------Method 1. Chain Ladder Method--------------#
library(ChainLadder)
# Rows: accident semesters (2018_S1 to 2024_S2)
# Columns: development periods (6, 12, 18, ..., 84 months)
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
  nrow = 14,
  byrow = TRUE
)
# Convert your matrix to a triangle object
triangle <- as.triangle(triangle_data)
# Print to check
print(triangle)
# ------------------------------------------------------------
# Traditional Chain Ladder Method
# ------------------------------------------------------------
# Add row and column names for clearer output
rownames(triangle) <- c(
  "2018_S1", "2018_S2",
  "2019_S1", "2019_S2",
  "2020_S1", "2020_S2",
  "2021_S1", "2021_S2",
  "2022_S1", "2022_S2",
  "2023_S1", "2023_S2",
  "2024_S1", "2024_S2"
)
colnames(triangle) <- paste0(seq(6, 84, by = 6), "m")
# Apply traditional deterministic Chain Ladder method
cl_model <- chainladder(triangle)
# Development factors
development_factors <- cl_model$f
# Completed cumulative triangle
full_triangle <- predict(cl_model)
cat("\n--- Full Cumulative Triangle ---\n")
print(full_triangle)
# Latest observed cumulative claims
latest <- getLatestCumulative(triangle)
# Ultimate claims estimates:
# correct extraction = last column of completed triangle
ultimate <- full_triangle[, ncol(full_triangle)]
# IBNR reserves
ibnr <- ultimate - latest
# Remove immaterial numerical noise
ibnr[abs(ibnr) < 1e-6] <- 0
# Results by accident period
results_by_origin <- data.frame(
  Accident_Period = rownames(triangle),
  Latest = as.numeric(latest),
  Ultimate = as.numeric(ultimate),
  IBNR = as.numeric(ibnr),
  row.names = NULL
)

# Totals, original scale
totals <- data.frame(
  Latest = sum(results_by_origin$Latest),
  Ultimate = sum(results_by_origin$Ultimate),
  IBNR = sum(results_by_origin$IBNR)
)

cat("\n--- Results by Accident Period, Original Scale ---\n")
print(results_by_origin)
cat("\n--- Totals, Original Scale ---\n")
print(totals)
# ------------------------------------------------------------
# Results in thousand ALL 
# ------------------------------------------------------------
results_by_origin_000 <- data.frame(
  Accident_Period = results_by_origin$Accident_Period,
  Latest = round(results_by_origin$Latest / 1000, 0),
  Ultimate = round(results_by_origin$Ultimate / 1000, 0),
  IBNR = round(results_by_origin$IBNR / 1000, 0),
  row.names = NULL
)

# For total row, calculate IBNR after rounding Latest and Ultimate
# Latest + IBNR = Ultimate
totals_000 <- data.frame(
  Latest = round(totals$Latest / 1000, 0),
  Ultimate = round(totals$Ultimate / 1000, 0)
)

totals_000$IBNR <- totals_000$Ultimate - totals_000$Latest
cat("\n--- Results by Accident Period, Thousand ALL ---\n")
print(results_by_origin_000)
cat("\n--- Totals, Thousand ALL ---\n")
print(totals_000)

# ============================================================
# Method 2: BORNHUETTER-FERGUSON METHOD
# ============================================================
# Description:
# This section applies the deterministic Bornhuetter-Ferguson
# method using the cumulative paid claims triangle, earned
# premiums and an a priori Expected Loss Ratio (ELR).
premiums <- c(
  136572286, 86234868,
  190102575, 152864071,
  168780826, 210152671,
  238006763, 186194472,
  286870652, 188831346,
  254656859, 253079992,
  329248107, 175794083
)

names(premiums) <- rownames(triangle)
# Expected Loss Ratio
ELR <- 0.75
# Expected ultimate claims
expected_ultimate <- premiums * ELR
# Convert triangle to matrix
triangle_matrix <- as.matrix(triangle)
#volume-weighted Chain Ladder development factors
n_dev <- ncol(triangle_matrix)
development_factors <- numeric (n_dev - 1)
for (j in 1:(n_dev - 1)) {
  valid_rows <- !is.na(triangle_matrix[, j]) & !is.na(triangle_matrix[, j + 1])
  development_factors[j] <- sum(triangle_matrix[valid_rows, j + 1]) /
    sum(triangle_matrix[valid_rows, j])
}
# Number of observed development periods by accident period
n_observed <- apply(triangle_matrix, 1, function(x) sum(!is.na(x)))
# CDF from latest observed development period to ultimate
cdf_to_ultimate <- sapply(n_observed, function(k) {
  
  if (k >= n_dev) {
    return(1)
  } else {
    return(prod(development_factors[k:(n_dev - 1)]))
  }
  
})

# Percentage reported and unreported proportion
percent_reported <- 1 / cdf_to_ultimate
unreported_proportion <- 1 - percent_reported

# Latest observed cumulative claims
latest <- getLatestCumulative(triangle)

# Bornhuetter-Ferguson reserve and ultimate
BF_IBNR <- expected_ultimate * unreported_proportion

BF_IBNR[abs(BF_IBNR) < 1e-6] <- 0

BF_ultimate <- latest + BF_IBNR

# Results by accident period
BF_results <- data.frame(
  Accident_Period = rownames(triangle),
  Latest = as.numeric(latest),
  Premiums = as.numeric(premiums),
  Expected_Ultimate = as.numeric(expected_ultimate),
  Unreported_Proportion = as.numeric(unreported_proportion),
  BF_IBNR = as.numeric(BF_IBNR),
  BF_Ultimate = as.numeric(BF_ultimate),
  row.names = NULL
)
# Totals
BF_totals <- data.frame(
  Latest = sum(BF_results$Latest),
  Premiums = sum(BF_results$Premiums),
  Expected_Ultimate = sum(BF_results$Expected_Ultimate),
  BF_IBNR = sum(BF_results$BF_IBNR),
  BF_Ultimate = sum(BF_results$BF_Ultimate)
)
# Results in thousand ALL
BF_results_000 <- data.frame(
  Accident_Period = BF_results$Accident_Period,
  Latest = round(BF_results$Latest / 1000, 0),
  Premiums = round(BF_results$Premiums / 1000, 0),
  Unreported_Proportion = round(BF_results$Unreported_Proportion, 5),
  Expected_Ultimate = round(BF_results$Expected_Ultimate / 1000, 0),
  BF_IBNR = round(BF_results$BF_IBNR / 1000, 0),
  BF_Ultimate = round(BF_results$BF_Ultimate / 1000, 0),
  row.names = NULL
)
BF_totals_000 <- data.frame(
  Latest = round(BF_totals$Latest / 1000, 0),
  Premiums = round(BF_totals$Premiums / 1000, 0),
  Expected_Ultimate = round(BF_totals$Expected_Ultimate / 1000, 0),
  BF_Ultimate = round(BF_totals$BF_Ultimate / 1000, 0)
)
BF_totals_000$BF_IBNR <- BF_totals_000$BF_Ultimate - BF_totals_000$Latest
# Print results
cat("\n--- Bornhuetter-Ferguson Results ---\n")
print(BF_results)
cat("\n--- Bornhuetter-Ferguson Totals ---\n")
print(BF_totals)

cat("\n--- Bornhuetter-Ferguson Results, in thousand ALL ---\n")
print(BF_results_000)

cat("\n--- Bornhuetter-Ferguson Totals, in thousand ALL ---\n")
print(BF_totals_000)

################################################################################
# Method 3.Double Chain Ladder (DCL) reserving method
# Reserves by origin period
# method based on two incremental triangles:
#   1. Incremental paid claims triangle
#   2. Incremental claim counts triangle
options(digits = 12)
#install.packages("DCL")
library(DCL)
library(dplyr)
library(tidyr)
library(ChainLadder)
library(clmplus)
################################################################################
# 1. Input data
# Semi-annual origin periods
origin_periods <- paste0(rep(2018:2024, each = 2), "_", rep(1:2, 7))
# Development periods
development_periods <- paste0("d", 0:13)

################################################################################
# 1.1 Incremental paid claims triangle
################################################################################

paid_triangle_inc <- matrix(
  c(
    69672343, 22708727,  157822,  93020,     0,      0,      0,  0,     0,     0,     0,     0,     0,     0,
    85370339, 21795046,  146568,      0,     0,      0,      0,  0,     0,     0,     0,     0,     0,    NA,
    99765261, 27313117,  340672,      0,     0,      0,      0,  0,     0,     0,     0,     0,    NA,    NA,
    87220584, 38945465,  430916,      0, 15520,      0,      0,  0,     0,     0,     0,    NA,    NA,    NA,
    75609735, 18383285,  217289,  22140,     0,      0,      0,  0,     0,     0,    NA,    NA,    NA,    NA,
    113677021, 38923559,  566224,    400,     0,      0,      0,  0,     0,    NA,    NA,    NA,    NA,    NA,
    144455980, 45354259,  756848, 200491,     0, 874976, 224651,  0,    NA,    NA,    NA,    NA,    NA,    NA,
    114860146, 36835581,  129222, 183926,  4000,  36784,      0, NA,    NA,    NA,    NA,    NA,    NA,    NA,
    145025581, 28817423, 1083627, 134420, 38500,      0,     NA, NA,    NA,    NA,    NA,    NA,    NA,    NA,
    128691752, 34997575,  135466,  86835,     0,     NA,     NA, NA,    NA,    NA,    NA,    NA,    NA,    NA,
    141608514, 27857928,   64289,   6400,    NA,     NA,     NA, NA,    NA,    NA,    NA,    NA,    NA,    NA,
    154641559, 30693999,  348480,     NA,    NA,     NA,     NA, NA,    NA,    NA,    NA,    NA,    NA,    NA,
    178455867, 30036188,      NA,     NA,    NA,     NA,     NA, NA,    NA,    NA,    NA,    NA,    NA,    NA,
    174183278,       NA,      NA,     NA,    NA,     NA,     NA, NA,    NA,    NA,    NA,    NA,    NA,    NA
  ),
  nrow = 14,
  byrow = TRUE,
  dimnames = list(origin_periods, development_periods)
)

################################################################################
# 1.2 Incremental reported claim counts triangle
################################################################################

count_triangle_inc <- matrix(
  c(
    5636, 162,  2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    5208, 283, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA,
    7698, 492, 34, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA,
    6247, 396, 37, 0, 1, 0, 0, 0, 0, 0, 0, NA, NA, NA,
    5429, 357, 19, 2, 0, 0, 0, 0, 0, 0, NA, NA, NA, NA,
    8029, 322, 24, 1, 0, 0, 0, 0, 0, NA, NA, NA, NA, NA,
    8841, 611, 17, 4, 0, 0, 0, 0, NA, NA, NA, NA, NA, NA,
    7991, 567,  6, 8, 1, 0, 0, NA, NA, NA, NA, NA, NA, NA,
    8714, 442, 40,11, 2, 0, NA, NA, NA, NA, NA, NA, NA, NA,
    8022, 749, 18, 1, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    8285, 861,  1, 1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    8900, 636,  0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    10200, 302, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    8936,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  nrow = 14,
  byrow = TRUE,
  dimnames = list(origin_periods, development_periods)
)

################################################################################
# 2. Data validation
################################################################################

validate_triangles <- function(paid_triangle, count_triangle) {
  
  stopifnot(is.matrix(paid_triangle))
  stopifnot(is.matrix(count_triangle))
  
  stopifnot(identical(dim(paid_triangle), dim(count_triangle)))
  stopifnot(identical(rownames(paid_triangle), rownames(count_triangle)))
  stopifnot(identical(colnames(paid_triangle), colnames(count_triangle)))
  
  stopifnot(all(paid_triangle[!is.na(paid_triangle)] >= 0))
  stopifnot(all(count_triangle[!is.na(count_triangle)] >= 0))
  
  invisible(TRUE)
}

validate_triangles(
  paid_triangle  = paid_triangle_inc,
  count_triangle = count_triangle_inc
)
################################################################################
# 3. Helper functions
################################################################################
remove_total_element <- function(x, expected_length) {
  
  if (is.null(x)) {
    stop("The prediction object does not contain the requested reserve vector.")
  }
  
  x <- as.numeric(x)
  
  if (length(x) > expected_length) {
    x <- head(x, expected_length)
  }
  
  return(x)
}

triangle_paid_to_date <- function(paid_triangle) {
  sum(rowSums(paid_triangle, na.rm = TRUE), na.rm = TRUE)
}
################################################################################
# 4. DCL model estimation
################################################################################

dcl_fit <- dcl.estimation(
  Xtriangle = paid_triangle_inc,
  Ntriangle = count_triangle_inc
)
################################################################################
# 5. DCL prediction by origin period
################################################################################
dcl_prediction_origin <- dcl.predict(
  dcl.par   = dcl_fit,
  Ntriangle = count_triangle_inc,
  Model     = 2,
  Tail      = TRUE,
  Tables    = FALSE,
  summ.by   = "row"
)
################################################################################
# 6. Extract reserve components by origin period
################################################################################
rbns_origin <- remove_total_element(
  x = dcl_prediction_origin$Rrbns,
  expected_length = length(origin_periods)
)
ibnr_origin <- remove_total_element(
  x = dcl_prediction_origin$Ribnr,
  expected_length = length(origin_periods)
)
total_origin <- remove_total_element(
  x = dcl_prediction_origin$Rtotal,
  expected_length = length(origin_periods)
)
################################################################################
# 7. Final reserve table by origin period
################################################################################
reserve_results_origin <- data.frame(
  Origin_Period  = origin_periods,
  RBNS           = rbns_origin,
  IBNR           = ibnr_origin,
  Total_Reserve  = total_origin,
  row.names = NULL
)
reserve_totals <- data.frame(
  Origin_Period  = "Total",
  RBNS           = sum(reserve_results_origin$RBNS, na.rm = TRUE),
  IBNR           = sum(reserve_results_origin$IBNR, na.rm = TRUE),
  Total_Reserve  = sum(reserve_results_origin$Total_Reserve, na.rm = TRUE),
  row.names = NULL
)
reserve_results_origin <- rbind(
  reserve_results_origin,
  reserve_totals
)
################################################################################
# 8. Paid-to-date and ultimate estimate
################################################################################
paid_to_date <- triangle_paid_to_date(paid_triangle_inc)
total_reserve <- reserve_totals$Total_Reserve
ultimate_estimate <- paid_to_date + total_reserve
summary_results <- data.frame(
  Measure = c(
    "Paid-to-date",
    "Total reserve",
    "Ultimate estimate"
  ),
  Amount = c(
    paid_to_date,
    total_reserve,
    ultimate_estimate
  )
)
################################################################################
# 9. Print results

cat("Double Chain Ladder reserves by origin period\n")
print(reserve_results_origin, row.names = FALSE)
cat("Summary results\n")
print(summary_results, row.names = FALSE)


#########################################################################
# Method 4. BORNHUETTER–DOUBLE CHAIN LADDER METHOD (BDCL)PREMIUM-BASED
suppressPackageStartupMessages({
  library(DCL)
})

if (exists("Xtriangle_inc") && !exists("paid_triangle_inc")) {
  paid_triangle_inc <- Xtriangle_inc
}

if (exists("counts_triangle") && !exists("count_triangle_inc")) {
  count_triangle_inc <- counts_triangle
}

premiums <- c(
  136572286, 86234868,
  190102575, 152864071,
  168780826, 210152671,
  238006763, 186194472,
  286870652, 188831346,
  254656859, 253079992,
  329248107, 175794083
)

# Expected Loss Ratio used as Bornhuetter–Ferguson prior.
ELR_prior <- 0.75
################################################################################
# 2. Input validation
validate_inputs <- function(paid_triangle_inc,
                            count_triangle_inc,
                            premiums,
                            ELR_prior) {
  
  if (!is.matrix(paid_triangle_inc)) {
    stop("paid_triangle_inc must be a matrix.")
  }
  
  if (!is.matrix(count_triangle_inc)) {
    stop("count_triangle_inc must be a matrix.")
  }
  
  if (!identical(dim(paid_triangle_inc), dim(count_triangle_inc))) {
    stop("The paid and count triangles must have the same dimensions.")
  }
  
  if (!identical(rownames(paid_triangle_inc), rownames(count_triangle_inc))) {
    stop("The paid and count triangles must have identical origin period names.")
  }
  
  if (!identical(colnames(paid_triangle_inc), colnames(count_triangle_inc))) {
    stop("The paid and count triangles must have identical development period names.")
  }
  
  if (length(premiums) != nrow(paid_triangle_inc)) {
    stop("The premiums vector must have the same length as the number of origin periods.")
  }
  
  if (any(is.na(premiums)) || any(premiums <= 0)) {
    stop("Premiums must be positive and non-missing.")
  }
  
  if (!is.numeric(ELR_prior) || length(ELR_prior) != 1 ||
      is.na(ELR_prior) || ELR_prior <= 0) {
    stop("ELR_prior must be a single positive numeric value.")
  }
  
  if (any(paid_triangle_inc[!is.na(paid_triangle_inc)] < 0)) {
    stop("paid_triangle_inc must not contain negative values.")
  }
  
  if (any(count_triangle_inc[!is.na(count_triangle_inc)] < 0)) {
    stop("count_triangle_inc must not contain negative values.")
  }
  
  invisible(TRUE)
}

validate_inputs(
  paid_triangle_inc  = paid_triangle_inc,
  count_triangle_inc = count_triangle_inc,
  premiums           = premiums,
  ELR_prior          = ELR_prior
)
################################################################################
# 3. Helper functions
# Convert an incremental triangle into a cumulative triangle.
incremental_to_cumulative <- function(triangle_inc) {
  triangle_inc <- as.matrix(triangle_inc)
  
  triangle_cum <- matrix(
    NA_real_,
    nrow = nrow(triangle_inc),
    ncol = ncol(triangle_inc),
    dimnames = dimnames(triangle_inc)
  )
  
  for (i in seq_len(nrow(triangle_inc))) {
    
    observed_cells <- which(!is.na(triangle_inc[i, ]))
    
    if (length(observed_cells) > 0) {
      triangle_cum[i, observed_cells] <- cumsum(triangle_inc[i, observed_cells])
    }
  }
  
  return(triangle_cum)
}
# Extract the latest observed cumulative value by origin period.

latest_cumulative_by_origin <- function(triangle_cum) {
  
  latest_values <- apply(triangle_cum, 1, function(x) {
    
    observed_values <- x[!is.na(x)]
    
    if (length(observed_values) == 0) {
      return(NA_real_)
    }
    
    return(tail(observed_values, 1))
  })
  
  return(as.numeric(latest_values))
}
# Calculate volume-weighted Chain Ladder age-to-age factors.
chain_ladder_factors <- function(triangle_cum) {
  
  triangle_cum <- as.matrix(triangle_cum)
  n_dev <- ncol(triangle_cum)
  
  factors <- rep(NA_real_, n_dev - 1)
  
  for (j in seq_len(n_dev - 1)) {
    
    valid_rows <- !is.na(triangle_cum[, j]) &
      !is.na(triangle_cum[, j + 1]) &
      triangle_cum[, j] > 0
    
    numerator <- sum(triangle_cum[valid_rows, j + 1], na.rm = TRUE)
    denominator <- sum(triangle_cum[valid_rows, j], na.rm = TRUE)
    
    if (denominator > 0) {
      factors[j] <- numerator / denominator
    } else {
      factors[j] <- 1
    }
  }
  
  return(factors)
}
# Calculate cumulative development factors from each development age to ultimate.

development_to_ultimate_factors <- function(age_to_age_factors) {
  
  n_dev <- length(age_to_age_factors) + 1
  
  cdf_to_ultimate <- rep(NA_real_, n_dev)
  cdf_to_ultimate[n_dev] <- 1
  
  for (j in (n_dev - 1):1) {
    
    factor_j <- age_to_age_factors[j]
    
    if (is.na(factor_j) || factor_j <= 0) {
      factor_j <- 1
    }
    
    cdf_to_ultimate[j] <- factor_j * cdf_to_ultimate[j + 1]
  }
  
  return(cdf_to_ultimate)
}
# Identify the latest observed development age by origin period.

latest_development_age <- function(triangle_cum) {
  
  latest_age <- apply(triangle_cum, 1, function(x) {
    
    observed_cells <- which(!is.na(x))
    
    if (length(observed_cells) == 0) {
      return(NA_integer_)
    }
    
    return(max(observed_cells))
  })
  
  return(as.integer(latest_age))
}
#extract a total value from a DCL prediction component.

extract_dcl_total <- function(x, component_name) {
  
  if (is.null(x)) {
    stop(paste("The DCL prediction object does not contain", component_name))
  }
  
  x <- as.numeric(x)
  
  if (length(x) == 0) {
    stop(paste(component_name, "is empty."))
  }
  
  return(as.numeric(tail(x, 1)))
}
################################################################################
# 4. Prepare cumulative triangles
################################################################################
paid_triangle_cum <- incremental_to_cumulative(paid_triangle_inc)
count_triangle_cum <- incremental_to_cumulative(count_triangle_inc)

origin_periods <- rownames(paid_triangle_inc)
n_origin <- nrow(paid_triangle_inc)

paid_to_date_origin <- latest_cumulative_by_origin(paid_triangle_cum)
################################################################################
# 5. Double Chain Ladder method
################################################################################

DCL_fit <- dcl.estimation(
  Xtriangle = paid_triangle_inc,
  Ntriangle = count_triangle_inc
)
DCL_prediction <- dcl.predict(
  dcl.par   = DCL_fit,
  Ntriangle = count_triangle_inc,
  Model     = 2,
  Tail      = TRUE,
  Tables    = FALSE
)

# DCL reserve components.
RBNS_DCL_total <- extract_dcl_total(
  x = DCL_prediction$Drbns,
  component_name = "Drbns"
)
IBNR_DCL_total <- extract_dcl_total(
  x = DCL_prediction$Dibnr,
  component_name = "Dibnr"
)
Total_DCL_total <- extract_dcl_total(
  x = DCL_prediction$Dtotal,
  component_name = "Dtotal"
)
################################################################################
# 6. Reported proportion from claim counts
################################################################################
count_age_to_age_factors <- chain_ladder_factors(count_triangle_cum)

count_cdf_to_ultimate <- development_to_ultimate_factors(
  age_to_age_factors = count_age_to_age_factors
)
count_latest_age <- latest_development_age(count_triangle_cum)

reported_proportion <- rep(NA_real_, n_origin)

for (i in seq_len(n_origin)) {
  
  age_i <- count_latest_age[i]
  
  if (is.na(age_i)) {
    reported_proportion[i] <- NA_real_
  } else {
    
    cdf_i <- count_cdf_to_ultimate[age_i]
    
    if (is.na(cdf_i) || cdf_i <= 0) {
      cdf_i <- 1
    }
    
    reported_proportion[i] <- min(1, 1 / cdf_i)
  }
}

unreported_proportion <- pmax(0, 1 - reported_proportion)
################################################################################
# 7. Bornhuetter–Ferguson prior based on premiums
################################################################################

Expected_Ultimate_BF <- premiums * ELR_prior

IBNR_BF_origin <- Expected_Ultimate_BF * unreported_proportion

IBNR_BF_origin[abs(IBNR_BF_origin) < 1e-6] <- 0

IBNR_BF_total <- sum(IBNR_BF_origin, na.rm = TRUE)
################################################################################
# 8. Allocation of DCL RBNS by origin period
# the total DCL RBNS component is allocated by origin period using the
# cumulative reported claim counts observed to date.

reported_counts_to_date <- latest_cumulative_by_origin(count_triangle_cum)
allocation_weights <- reported_counts_to_date
allocation_weights[is.na(allocation_weights)] <- 0
if (sum(allocation_weights, na.rm = TRUE) <= 0) {
  allocation_weights <- paid_to_date_origin
  allocation_weights[is.na(allocation_weights)] <- 0
}

allocation_weights <- allocation_weights / sum(allocation_weights, na.rm = TRUE)

RBNS_DCL_origin <- RBNS_DCL_total * allocation_weights

################################################################################
# 9. BDCL reserve and ultimate by origin period
################################################################################
Total_Reserve_BDCL_P_origin <- RBNS_DCL_origin + IBNR_BF_origin

Ultimate_BDCL_P_origin <- paid_to_date_origin + Total_Reserve_BDCL_P_origin

################################################################################
# 10. Final result table by origin period
################################################################################

BDCL_P_results_origin <- data.frame(
  Origin_Period            = origin_periods,
  Premium                  = as.numeric(premiums),
  ELR_Prior                = rep(ELR_prior, n_origin),
  Paid_to_Date             = as.numeric(paid_to_date_origin),
  Reported_Proportion      = as.numeric(reported_proportion),
  Unreported_Proportion    = as.numeric(unreported_proportion),
  Expected_Ultimate_BF     = as.numeric(Expected_Ultimate_BF),
  RBNS_DCL                 = as.numeric(RBNS_DCL_origin),
  IBNR_BF                  = as.numeric(IBNR_BF_origin),
  Total_Reserve_BDCL_P     = as.numeric(Total_Reserve_BDCL_P_origin),
  Ultimate_BDCL_P          = as.numeric(Ultimate_BDCL_P_origin),
  row.names = NULL
)

BDCL_P_total_row <- data.frame(
  Origin_Period            = "Total",
  Premium                  = sum(BDCL_P_results_origin$Premium, na.rm = TRUE),
  ELR_Prior                = ELR_prior,
  Paid_to_Date             = sum(BDCL_P_results_origin$Paid_to_Date, na.rm = TRUE),
  Reported_Proportion      = NA_real_,
  Unreported_Proportion    = NA_real_,
  Expected_Ultimate_BF     = sum(BDCL_P_results_origin$Expected_Ultimate_BF, na.rm = TRUE),
  RBNS_DCL                 = sum(BDCL_P_results_origin$RBNS_DCL, na.rm = TRUE),
  IBNR_BF                  = sum(BDCL_P_results_origin$IBNR_BF, na.rm = TRUE),
  Total_Reserve_BDCL_P     = sum(BDCL_P_results_origin$Total_Reserve_BDCL_P, na.rm = TRUE),
  Ultimate_BDCL_P          = sum(BDCL_P_results_origin$Ultimate_BDCL_P, na.rm = TRUE),
  row.names = NULL
)

BDCL_P_results_origin <- rbind(
  BDCL_P_results_origin,
  BDCL_P_total_row
)
################################################################################
# 11. Results in thousand monetary units
################################################################################

BDCL_P_results_origin_000 <- BDCL_P_results_origin

amount_columns <- c(
  "Premium",
  "Paid_to_Date",
  "Expected_Ultimate_BF",
  "RBNS_DCL",
  "IBNR_BF",
  "Total_Reserve_BDCL_P",
  "Ultimate_BDCL_P"
)

BDCL_P_results_origin_000[, amount_columns] <- round(
  BDCL_P_results_origin_000[, amount_columns] / 1000,
  0
)

BDCL_P_results_origin_000$ELR_Prior <- round(
  BDCL_P_results_origin_000$ELR_Prior,
  6
)

BDCL_P_results_origin_000$Reported_Proportion <- round(
  BDCL_P_results_origin_000$Reported_Proportion,
  6
)

BDCL_P_results_origin_000$Unreported_Proportion <- round(
  BDCL_P_results_origin_000$Unreported_Proportion,
  6
)
################################################################################
# 12. Summary tables
BDCL_P_summary <- data.frame(
  Measure = c(
    "Paid-to-date",
    "Expected ultimate BF",
    "RBNS reserve from DCL",
    "IBNR reserve from BF",
    "Total BDCL-P reserve",
    "Ultimate BDCL-P estimate"
  ),
  Amount = c(
    sum(paid_to_date_origin, na.rm = TRUE),
    sum(Expected_Ultimate_BF, na.rm = TRUE),
    sum(RBNS_DCL_origin, na.rm = TRUE),
    sum(IBNR_BF_origin, na.rm = TRUE),
    sum(Total_Reserve_BDCL_P_origin, na.rm = TRUE),
    sum(Ultimate_BDCL_P_origin, na.rm = TRUE)
  ),
  row.names = NULL
)

BDCL_P_summary_000 <- BDCL_P_summary
BDCL_P_summary_000$Amount <- round(BDCL_P_summary_000$Amount / 1000, 0)
################################################################################
# 13. Consistency checks

stopifnot(
  abs(sum(RBNS_DCL_origin, na.rm = TRUE) - RBNS_DCL_total) < 1e-6
)

stopifnot(
  abs(sum(IBNR_BF_origin, na.rm = TRUE) - IBNR_BF_total) < 1e-6
)

stopifnot(
  abs(
    sum(Total_Reserve_BDCL_P_origin, na.rm = TRUE) -
      (
        sum(RBNS_DCL_origin, na.rm = TRUE) +
          sum(IBNR_BF_origin, na.rm = TRUE)
      )
  ) < 1e-6
)

stopifnot(
  abs(
    sum(Ultimate_BDCL_P_origin, na.rm = TRUE) -
      (
        sum(paid_to_date_origin, na.rm = TRUE) +
          sum(Total_Reserve_BDCL_P_origin, na.rm = TRUE)
      )
  ) < 1e-6
)

################################################################################
# 14. Print results
cat("Premium-based Bornhuetter–Double Chain Ladder results by origin period\n")
print(BDCL_P_results_origin, row.names = FALSE)
cat("Premium-based Bornhuetter–Double Chain Ladder results by origin period")
cat(", amounts in thousands\n")
print(BDCL_P_results_origin_000, row.names = FALSE)
cat("Premium-based Bornhuetter–Double Chain Ladder summary results\n")
print(BDCL_P_summary, row.names = FALSE)
cat("Premium-based Bornhuetter–Double Chain Ladder summary results")
cat(", amounts in thousands\n")
print(BDCL_P_summary_000, row.names = FALSE)



############################################################
###              STOCHASTIC METHODS
############################################################
# This section includes stochastic reserving methods,
# where, in addition to estimating reserves, the uncertainty
# of the prediction is also measured.

##### Method 5.Mack Chain-Ladder method ##
# The object 'triangle' is the cumulative paid claims triangle
# Run Mack Chain-Ladder with Mack estimation of standard errors
mack_result <- MackChainLadder(
  Triangle = triangle,
  est.sigma = "Mack"
)
# Summary of Mack results:
# Latest, development-to-date, ultimate claims, IBNR,
# Mack standard error and coefficient of variation
mack_summary <- summary(mack_result)
print(mack_summary)

# Diagnostic plots
plot(mack_result)
# Development factors
mack_development_factors <- mack_result$f
print(mack_development_factors)

# Cumulative development factors
mack_cumulative_factors <- rev(cumprod(rev(mack_result$f)))
cat("\nCumulative Development Factors:\n")
print(mack_cumulative_factors)

# Full triangle with Mack forecasts
mack_full_triangle <- mack_result$FullTriangle
print(mack_full_triangle)

# Mack standard errors of IBNR
mack_standard_errors <- mack_result$Mack.S.E
print(mack_standard_errors)

############## Method 6. Munich Chain Ladder Method #####
library(ChainLadder)
# ------------------------------------------------------------
devs    <- paste0("Dev", seq(6, 84, 6))
origins <- c("2018_1","2018_2","2019_1","2019_2","2020_1","2020_2",
             "2021_1","2021_2","2022_1","2022_2","2023_1","2023_2",
             "2024_1","2024_2")

# ---- INCURRED (cumulative) ----
incurred <- matrix(
  byrow = TRUE, nrow = length(origins),
  c(
    102700465,126789437,127130337,127223357,127223357,127223357,127223357,127223357,127223357,127223357,127223357,127223357,127223357,127223357,
    118630052,143874164,144094016,144094016,144094016,144094016,144094016,144094016,144094016,144094016,144094016,144094016,144094016,NA,
    147931805,177099629,177607937,177607937,177607937,177607937,177607937,177607937,177607937,177607937,177607937,177607937,NA,NA,
    134951237,179045040,179638392,179638392,179669432,179669432,179669432,179669432,179669432,179669432,179669432,NA,NA,NA,
    109719349,129208874,129523717,129550261,129550261,129550261,129550261,129550261,129550261,129550261,NA,NA,NA,NA,
    172595154,214137036,214826643,214827043,214827043,214827043,214827043,214827043,214827043,NA,NA,NA,NA,NA,
    218652356,272431302,275659815,277927894,279806694,280681670,280906321,280906321,NA,NA,NA,NA,NA,NA,
    178396460,225182241,225697796,225886246,225894246,225931030,225931030,NA,NA,NA,NA,NA,NA,NA,
    210154517,242959387,245418238,245704050,245745050,245745050,NA,NA,NA,NA,NA,NA,NA,NA,
    173818146,212391526,212761028,213054432,213054432,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    191682554,222744711,223630496,223636896,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    209212481,245296561,245645041,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    253139091,290922272,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    213451798,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
  )
)
rownames(incurred) <- origins
colnames(incurred) <- devs
# ---- PAID (cumulative) ----
paid <- matrix(
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
  nrow = 14, byrow = TRUE
)
rownames(paid) <- origins
colnames(paid) <- devs
# ------------------------------------------------------------
# Convert matrices to triangle objects
PAID <- as.triangle(paid)
INC  <- as.triangle(incurred)

# Munich Chain Ladder
MCL <- MunichChainLadder(
  Paid = PAID,
  Incurred = INC,
  est.sigmaP = "Mack",
  est.sigmaI = "Mack"
)

# Results by origin period
MCLPaid_mat     <- as.matrix(MCL$MCLPaid)
MCLIncurred_mat <- as.matrix(MCL$MCLIncurred)

Latest_Paid     <- as.numeric(getLatestCumulative(PAID))
Latest_Incurred <- as.numeric(getLatestCumulative(INC))

Ultimate_Paid     <- as.numeric(MCLPaid_mat[, ncol(MCLPaid_mat)])
Ultimate_Incurred <- as.numeric(MCLIncurred_mat[, ncol(MCLIncurred_mat)])

Reserve_Paid     <- Ultimate_Paid - Latest_Paid
Reserve_Incurred <- Ultimate_Incurred - Latest_Incurred

# ------------------------------------------------------------
# Results by origin period - original values
# ------------------------------------------------------------

Munich_By_Period <- data.frame(
  Period = rownames(PAID),
  Paid_Cumulative = Latest_Paid,
  Incurred_Cumulative = Latest_Incurred,
  Current_PI_Ratio = Latest_Paid / Latest_Incurred,
  Ultimate_Paid = Ultimate_Paid,
  Ultimate_Incurred = Ultimate_Incurred,
  Final_PI_Ratio = Ultimate_Paid / Ultimate_Incurred,
  Reserve_Paid = Reserve_Paid,
  Reserve_Incurred = Reserve_Incurred
)

print(Munich_By_Period)

# ------------------------------------------------------------
# Results by origin period - rounded in thousand ALL
# ------------------------------------------------------------

Munich_By_Period_000 <- data.frame(
  Period = rownames(PAID),
  Paid_Cumulative = round(Latest_Paid / 1000, 0),
  Incurred_Cumulative = round(Latest_Incurred / 1000, 0),
  Current_PI_Ratio = round(Latest_Paid / Latest_Incurred, 3),
  Ultimate_Paid = round(Ultimate_Paid / 1000, 0),
  Ultimate_Incurred = round(Ultimate_Incurred / 1000, 0),
  Final_PI_Ratio = round(Ultimate_Paid / Ultimate_Incurred, 3),
  Reserve_Paid = round(Reserve_Paid / 1000, 0),
  Reserve_Incurred = round(Reserve_Incurred / 1000, 0)
)

print(Munich_By_Period_000)

Munich_Total <- data.frame(
  Measure = c("Ultimate", "Reserve"),
  Paid = c(
    sum(Ultimate_Paid, na.rm = TRUE),
    sum(Reserve_Paid, na.rm = TRUE)
  ),
  Incurred = c(
    sum(Ultimate_Incurred, na.rm = TRUE),
    sum(Reserve_Incurred, na.rm = TRUE)
  )
)

print(Munich_Total)
# ------------------------------------------------------------
# Total results - rounded in thousand ALL
# ------------------------------------------------------------
Total_Latest_Paid_000     <- round(sum(Latest_Paid, na.rm = TRUE) / 1000, 0)
Total_Latest_Incurred_000 <- round(sum(Latest_Incurred, na.rm = TRUE) / 1000, 0)

Total_Ultimate_Paid_000     <- round(sum(Ultimate_Paid, na.rm = TRUE) / 1000, 0)
Total_Ultimate_Incurred_000 <- round(sum(Ultimate_Incurred, na.rm = TRUE) / 1000, 0)

Total_Reserve_Paid_000     <- Total_Ultimate_Paid_000 - Total_Latest_Paid_000
Total_Reserve_Incurred_000 <- Total_Ultimate_Incurred_000 - Total_Latest_Incurred_000

Munich_Total_000 <- data.frame(
  Measure = c("Cumulative", "Ultimate", "Reserve"),
  Paid = c(
    Total_Latest_Paid_000,
    Total_Ultimate_Paid_000,
    Total_Reserve_Paid_000
  ),
  Incurred = c(
    Total_Latest_Incurred_000,
    Total_Ultimate_Incurred_000,
    Total_Reserve_Incurred_000
  )
)

print(Munich_Total_000)
# Diagnostic plots
plot(MCL)
# ------------------------------------------------------------
# Standard deviation and CV
# ------------------------------------------------------------
SD_Paid     <- as.numeric(MCL$MackPaid$Mack.S.E)
SD_Incurred <- as.numeric(MCL$MackIncurred$Mack.S.E)

CV_Paid     <- ifelse(Reserve_Paid == 0, NA, SD_Paid / Reserve_Paid)
CV_Incurred <- ifelse(Reserve_Incurred == 0, NA, SD_Incurred / Reserve_Incurred)

Munich_Total_SD_CV <- data.frame(
  Measure = c("Ultimate", "Reserve", "Standard Deviation", "CV"),
  Paid = c(
    round(sum(Ultimate_Paid, na.rm = TRUE) / 1000, 0),
    Total_Reserve_Paid_000,
    round(sqrt(sum(SD_Paid^2, na.rm = TRUE)) / 1000, 0),
    round(sqrt(sum(SD_Paid^2, na.rm = TRUE)) / sum(Reserve_Paid, na.rm = TRUE), 2)
  ),
  Incurred = c(
    round(sum(Ultimate_Incurred, na.rm = TRUE) / 1000, 0),
    Total_Reserve_Incurred_000,
    round(sqrt(sum(SD_Incurred^2, na.rm = TRUE)) / 1000, 0),
    round(sqrt(sum(SD_Incurred^2, na.rm = TRUE)) / sum(Reserve_Incurred, na.rm = TRUE), 2)
  )
)

print(Munich_Total_SD_CV)

#------------------------------------------------------------#
#    Clark's Methods      #
#--------------------------------------------
#---------------Method 7. Clark LDF Method with Loglogistic Distribution----------#

triangle <- as.triangle(triangle_data)

cldf <- ClarkLDF(
  Triangle = triangle,
  cumulative = TRUE,
  maxage = 84,
  adol = FALSE,
  adol.age = NULL,
  origin.width = NULL,
  G = "loglogistic"
)
# Model results
cldf
# Parameters of the loglogistic growth curve
cldf$THETAG
# Diagnostic plot
graphics.off()
par(
  mfrow = c(1, 1),
  mar = c(4, 4, 2, 1) + 0.1,
  oma = c(0, 0, 0, 0)
)

plot(cldf)


#---------------Method 8. Clark Cape Cod Method with Loglogistic Distribution----------#

premium_vector <- c(
  136572286, 86234868, 190102575, 152864071,
  168780826, 210152671, 238006763, 186194472,
  286870652, 188831346, 254656859, 253079992,
  329248107, 175794083
)

colnames(triangle_data) <- seq(6, 6 * ncol(triangle_data), by = 6)

rownames(triangle_data) <- c(
  "2018_1", "2018_2",
  "2019_1", "2019_2",
  "2020_1", "2020_2",
  "2021_1", "2021_2",
  "2022_1", "2022_2",
  "2023_1", "2023_2",
  "2024_1", "2024_2"
)

triangle <- as.triangle(triangle_data)

cccll <- ClarkCapeCod(
  Triangle = triangle,
  Premium = premium_vector,
  cumulative = TRUE,
  maxage = 84,
  adol = FALSE,
  adol.age = NULL,
  origin.width = NULL,
  G = "loglogistic"
)

summary(cccll)

plot(cccll)

#---------------------------------------------------
#          Bootstrap methods
#--------------------------------------------------
##------ Method 9. Bootstrap CL Method------------#########

# Reproducibility
set.seed(123)

# Triangle object
triangle <- as.triangle(triangle_data)

# Bootstrap simulation
boot_result <- BootChainLadder(
  triangle,
  R = 10000,
  process.distr = "gamma"
)

# Summary
summary(boot_result)

# Quantiles for VaR
VaR_probs <- c(0.75, 0.90, 0.95, 0.995)

VaR_values <- quantile(boot_result$IBNR.Totals, VaR_probs)
# Print
cat("Bootstrap IBNR VaR Quantiles:\n")
print(data.frame(
  Probability = VaR_probs,
  VaR = round(VaR_values, 2)
))

quantile(boot_result, c(0.75,0.95,0.99, 0.995))

# Plot ECDF
plot(
  ecdf(boot_result$IBNR.Totals),
  main = "ECDF of Bootstrap Total IBNR",
  xlab = "IBNR Reserve",
  ylab = "Empirical Probability",
  col = "blue",
  lwd = 2
)
abline(v = VaR_values, col = "red", lty = 2)
legend(
  "bottomright",
  legend = paste0("VaR ", VaR_probs * 100, "%"),
  col = "red",
  lty = 2,
  cex = 0.8
)

hist(
  boot_result$IBNR.Totals,
  breaks=30,
  col="gray",
  border="white",
  main="Histogram of Total IBNR",
  xlab="Total IBNR"
)
plot(
  ecdf(boot_result$IBNR.Totals),
  main="ECDF of Total IBNR",
  xlab="Total IBNR",
  ylab="F(x)"
)
grid()

plot(boot_result)
graphics.off()
par(mfrow = c(1,1), mar = c(5,4,2,1) + 0.1, oma = c(0,0,0,0))
plot(boot_result)

summary(boot_result$IBNR.Totals)

################################################################################
# Method 10. Bootstrap Double Chain Ladder Method (Bootstrap DCL)
################################################################################
options(digits = 12)

library(DCL)

################################################################################
# 1 Data validation
################################################################################

stopifnot(identical(dim(paid_triangle_inc), dim(count_triangle_inc)))
stopifnot(identical(rownames(paid_triangle_inc), rownames(count_triangle_inc)))
stopifnot(identical(colnames(paid_triangle_inc), colnames(count_triangle_inc)))

origin_labels <- rownames(paid_triangle_inc)

################################################################################
# 2 Fit the DCL model
################################################################################

dcl_fit <- dcl.estimation(
  Xtriangle = paid_triangle_inc,
  Ntriangle = count_triangle_inc
)

################################################################################
# 3 Classical DCL point estimates by origin period

dcl_cell_prediction <- dcl.predict(
  dcl.par   = dcl_fit,
  Ntriangle = count_triangle_inc,
  Model     = 2,
  Tail      = TRUE,
  Tables    = FALSE,
  summ.by   = "cell"
)

rbns_origin_point <- rowSums(dcl_cell_prediction$Xrbns, na.rm = TRUE)
ibnr_origin_point <- rowSums(dcl_cell_prediction$Xibnr, na.rm = TRUE)
total_origin_point <- rowSums(dcl_cell_prediction$Xtotal, na.rm = TRUE)

rbns_total_point <- sum(rbns_origin_point, na.rm = TRUE)
ibnr_total_point <- sum(ibnr_origin_point, na.rm = TRUE)
total_total_point <- sum(total_origin_point, na.rm = TRUE)

################################################################################
# 4 Run the DCL bootstrap
################################################################################

number_of_simulations <- 10000

set.seed(123)

dcl_bootstrap <- dcl.boot(
  dcl.par   = dcl_fit,
  Ntriangle = count_triangle_inc,
  boot.type = 2,
  B         = number_of_simulations,
  Tail      = TRUE,
  Tables    = FALSE
)

################################################################################
# 5 Extract bootstrap reserves by origin period
################################################################################
# The DCL bootstrap output is stored in three-dimensional arrays:
#   origin period x development/calendar dimension x bootstrap simulation

extract_origin_bootstrap <- function(bootstrap_array, origin_labels) {
  
  array_dimensions <- dim(bootstrap_array)
  
  if (length(array_dimensions) != 3) {
    stop("The bootstrap object must be a three-dimensional array.")
  }
  
  number_of_origins <- array_dimensions[1]
  number_of_bootstrap_runs <- array_dimensions[3]
  
  if (number_of_origins != length(origin_labels)) {
    stop("The number of origin labels does not match the bootstrap array.")
  }
  
  bootstrap_matrix <- matrix(
    NA_real_,
    nrow = number_of_bootstrap_runs,
    ncol = number_of_origins
  )
  
  for (b in seq_len(number_of_bootstrap_runs)) {
    bootstrap_matrix[b, ] <- rowSums(bootstrap_array[, , b], na.rm = TRUE)
  }
  
  colnames(bootstrap_matrix) <- origin_labels
  
  return(bootstrap_matrix)
}

rbns_boot_origin <- extract_origin_bootstrap(
  bootstrap_array = dcl_bootstrap$array.rbns.boot,
  origin_labels   = origin_labels
)

ibnr_boot_origin <- extract_origin_bootstrap(
  bootstrap_array = dcl_bootstrap$array.ibnr.boot,
  origin_labels   = origin_labels
)

total_boot_origin <- rbns_boot_origin + ibnr_boot_origin

################################################################################
# 6 Add total column across all origin periods
################################################################################

rbns_boot_origin <- cbind(
  rbns_boot_origin,
  Total = rowSums(rbns_boot_origin, na.rm = TRUE)
)

ibnr_boot_origin <- cbind(
  ibnr_boot_origin,
  Total = rowSums(ibnr_boot_origin, na.rm = TRUE)
)

total_boot_origin <- cbind(
  total_boot_origin,
  Total = rowSums(total_boot_origin, na.rm = TRUE)
)

################################################################################
# 7 Create bootstrap summary tables
################################################################################
# Percentiles used in the dissertation:
#   P1  = 1st percentile
#   P5  = 5th percentile
#   P50 = 50th percentile, median
#   P95 = 95th percentile
#   P99 = 99th percentile

make_bootstrap_summary <- function(bootstrap_matrix, point_estimates, origin_labels) {
  
  percentile_matrix <- t(apply(
    bootstrap_matrix,
    2,
    quantile,
    probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
    na.rm = TRUE,
    type = 7
  ))
  
  summary_table <- data.frame(
    Period = c(origin_labels, "Total"),
    Claim_Reserve = point_estimates,
    Simulated_Mean_Reserve = colMeans(bootstrap_matrix, na.rm = TRUE),
    Reserve_Standard_Deviation = apply(bootstrap_matrix, 2, sd, na.rm = TRUE),
    P1 = percentile_matrix[, 1],
    P5 = percentile_matrix[, 2],
    P50_Median = percentile_matrix[, 3],
    P95 = percentile_matrix[, 4],
    P99 = percentile_matrix[, 5],
    row.names = NULL,
    check.names = FALSE
  )
  
  return(summary_table)
}

rbns_bootstrap_summary <- make_bootstrap_summary(
  bootstrap_matrix = rbns_boot_origin,
  point_estimates  = c(rbns_origin_point, rbns_total_point),
  origin_labels    = origin_labels
)

ibnr_bootstrap_summary <- make_bootstrap_summary(
  bootstrap_matrix = ibnr_boot_origin,
  point_estimates  = c(ibnr_origin_point, ibnr_total_point),
  origin_labels    = origin_labels
)

total_bootstrap_summary <- make_bootstrap_summary(
  bootstrap_matrix = total_boot_origin,
  point_estimates  = c(total_origin_point, total_total_point),
  origin_labels    = origin_labels
)

################################################################################
# 8 Format monetary results in thousands
################################################################################
format_in_thousands <- function(summary_table) {
  
  formatted_table <- summary_table
  numeric_columns <- sapply(formatted_table, is.numeric)
  
  formatted_table[numeric_columns] <- lapply(
    formatted_table[numeric_columns],
    function(x) round(x / 1000)
  )
  
  return(formatted_table)
}

rbns_bootstrap_summary_thousand <- format_in_thousands(rbns_bootstrap_summary)
ibnr_bootstrap_summary_thousand <- format_in_thousands(ibnr_bootstrap_summary)
total_bootstrap_summary_thousand <- format_in_thousands(total_bootstrap_summary)

################################################################################
# 9 Print final results
cat("TOTAL RESERVE BOOTSTRAP SUMMARY BY ORIGIN PERIOD\n")
cat("Percentiles: P1, P5, P50, P95, P99\n")
print(total_bootstrap_summary_thousand, row.names = FALSE)

cat("\n=============================================================================\n")
cat("RBNS BOOTSTRAP SUMMARY BY ORIGIN PERIOD\n")
cat("Percentiles: P1, P5, P50, P95, P99\n")

print(rbns_bootstrap_summary_thousand, row.names = FALSE)
cat("Percentiles: P1, P5, P50, P95, P99\n")
print(ibnr_bootstrap_summary_thousand, row.names = FALSE)


########################################
paid_to_date_by_origin <- rowSums(paid_triangle_inc, na.rm = TRUE)
paid_to_date_total <- sum(paid_to_date_by_origin, na.rm = TRUE)

paid_to_date <- c(
  paid_to_date_by_origin,
  paid_to_date_total
)

total_bootstrap_summary$Paid_to_Date <- paid_to_date

total_bootstrap_summary$Ultimate_Value <- 
  total_bootstrap_summary$Paid_to_Date +
  total_bootstrap_summary$Claim_Reserve

total_bootstrap_summary$Reserve_Coefficient_of_Variation <- ifelse(
  total_bootstrap_summary$Claim_Reserve > 0,
  total_bootstrap_summary$Reserve_Standard_Deviation /
    total_bootstrap_summary$Claim_Reserve,
  NA_real_
)

total_bootstrap_summary <- total_bootstrap_summary[, c(
  "Period",
  "Paid_to_Date",
  "Ultimate_Value",
  "Claim_Reserve",
  "Simulated_Mean_Reserve",
  "Reserve_Standard_Deviation",
  "Reserve_Coefficient_of_Variation",
  "P1",
  "P5",
  "P50_Median",
  "P95",
  "P99"
)]

print(total_bootstrap_summary, row.names = FALSE)

###################  END ########################


