#' compare area and partial area based on logit transformation for ordinal data and unpaired design
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param e1 lower bound of FPR rate for partial area
#' @param e2 upper bound of FPR rate for partial area
#' @param alpha significance level
#' @return list containing test statistics and whether the areas or partial areas are significantly different using different transformations
#' example
#' data1_no<-c(38,25,15,19,4)
#' data1_yes<-c(1,2,3,14,42)
#' data2_no<-c(70,7,5,7,12)
#' data2_yes<-c(8,2,2,2,48)
#' area_logit_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0, 0.2, 0.05)
#' @export
area_logit_unpaired <- function(data1_no, data1_yes, data2_no, data2_yes, e1, e2, alpha) {
  estimation1 <- params_estimate_ordinal_unpaired(data1_no, data1_yes)
  estimation2 <- params_estimate_ordinal_unpaired(data2_no, data2_yes)

  # Full area
  a1 <- estimation1$a
  b1 <- estimation1$b
  var_a1 <- estimation1$cov_matrix[1, 1]
  var_b1 <- estimation1$cov_matrix[2, 2]
  covar_a1b1 <- estimation1$cov_matrix[1, 2]

  a2 <- estimation2$a
  b2 <- estimation2$b
  var_a2 <- estimation2$cov_matrix[1, 1]
  var_b2 <- estimation2$cov_matrix[2, 2]
  covar_a2b2 <- estimation2$cov_matrix[1, 2]

  A1 <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$area.full
  A2 <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$area.full

  varlogitA1 <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_full / (A1 * (1 - A1))^2
  varlogitA2 <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_full / (A2 * (1 - A2))^2

  Z_full <- (log(A1 / (1 - A1)) - log(A2 / (1 - A2))) / sqrt(varlogitA1 + varlogitA2)

  full_area_result <- ifelse(abs(Z_full) < qnorm(1 - alpha / 2),
                             "Using logit transformation, the full areas are not significantly different",
                             "Using logit transformation, the full areas are significantly different")

  # Partial area (Logit Transformation)
  A1_partial <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$area.partial
  A2_partial <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$area.partial

  varlogitA1_partial <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_partial / (A1_partial * (1 - A1_partial))^2
  varlogitA2_partial <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_partial / (A2_partial * (1 - A2_partial))^2

  Z_partial <- (log(A1_partial / (1 - A1_partial)) - log(A2_partial / (1 - A2_partial))) / sqrt(varlogitA1_partial + varlogitA2_partial)

  partial_area_logit_result <- ifelse(abs(Z_partial) < qnorm(1 - alpha / 2),
                                      "Using logit transformation, the partial areas are not significantly different",
                                      "Using logit transformation, the partial areas are significantly different")

  # Partial area (McClish Transformation)
  Amx <- e2 - e1
  varMcClishA1 <- 4 * Amx^2 * area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_partial / (Amx^2 - A1_partial^2)^2
  varMcClishA2 <- 4 * Amx^2 * area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_partial / (Amx^2 - A2_partial^2)^2

  Z_partial_McClish <- (log((Amx + A1_partial) / (Amx - A1_partial)) - log((Amx + A2_partial) / (Amx - A2_partial))) / sqrt(varMcClishA1 + varMcClishA2)

  partial_area_mcclish_result <- ifelse(abs(Z_partial_McClish) < qnorm(1 - alpha / 2),
                                        "Using McClish transformation, the partial areas are not significantly different",
                                        "Using McClish transformation, the partial areas are significantly different")

  # Return results
  result <- list(
    test_statistic_full = Z_full,
    full_area_result = full_area_result,
    test_statistic_partial_logit = Z_partial,
    partial_area_logit_result = partial_area_logit_result,
    test_statistic_partial_mcclish = Z_partial_McClish,
    partial_area_mcclish_result = partial_area_mcclish_result
  )

  return(result)
}





#' compare area and partial area based on logit transformation for ordinal data and paired design
#' @param data0 a matrix representing paired test results for undiseased subjects
#' @param data1 a matrix representing paired test results for diseased subjects
#' @param e1 lower bound of FPR rate for partial area
#' @param e2 upper bound of FPR rate for partial area
#' @param alpha significance level
#' @return list containing test statistics and whether the areas or partial areas are significantly different using different transformations
#' example
#' data0<-matrix(c(36,20,8,6,0,0,3,3,1,0,0,1,0,4,0,0,0,2,3,2,2,1,2,5,2),nrow=5,byrow=TRUE)
#' data1<-matrix(c(1,2,2,2,1,0,0,0,0,2,0,0,0,2,0,0,0,0,2,0,0,0,1,8,39),nrow=5,byrow=TRUE)
#' area_logit_paired(data0, data1, 0, 0.2, 0.05)
#' @export
area_logit_paired <- function(data0, data1, e1, e2, alpha) {
  estimation <- params_estimate_ordinal_paired(data0, data1)

  # Full area
  a1 <- estimation$a1
  b1 <- estimation$b1
  var_a1 <- estimation$cov_matrix[1, 1]
  var_b1 <- estimation$cov_matrix[2, 2]
  covar_a1b1 <- estimation$cov_matrix[1, 2]

  a2 <- estimation$a2
  b2 <- estimation$b2
  var_a2 <- estimation$cov_matrix[3, 3]
  var_b2 <- estimation$cov_matrix[4, 4]
  covar_a2b2 <- estimation$cov_matrix[3, 4]

  A1 <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$area.full
  A2 <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$area.full

  varlogitA1 <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_full / (A1 * (1 - A1))^2
  varlogitA2 <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_full / (A2 * (1 - A2))^2

  # Delta method for covariance matrix
  grad <- matrix(0, nrow = 2, ncol = 4)
  grad[1, 1] <- dnorm(a1 / sqrt(1 + b1^2)) / (sqrt(1 + b1^2))
  grad[1, 2] <- b1 * dnorm(a1 / sqrt(1 + b1^2)) / sqrt((1 + b1^2)^3)
  grad[2, 3] <- dnorm(a2 / sqrt(1 + b2^2)) / (sqrt(1 + b2^2))
  grad[2, 4] <- b2 * dnorm(a2 / sqrt(1 + b2^2)) / sqrt((1 + b2^2)^3)

  Delta_matrix <- grad %*% estimation$cov_matrix %*% t(grad)
  covlogit <- Delta_matrix[1, 2] / (A1 * A2 * (1 - A1) * (1 - A2))

  Z_full <- (log(A1 / (1 - A1)) - log(A2 / (1 - A2))) / sqrt(varlogitA1 + varlogitA2 - 2 * covlogit)

  full_area_result <- ifelse(abs(Z_full) < qnorm(1 - alpha / 2),
                             "Using logit transformation, the full areas are not significantly different",
                             "Using logit transformation, the full areas are significantly different")

  # Partial area (Logit Transformation)
  A1_partial <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$area.partial
  A2_partial <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$area.partial

  varlogitA1_partial <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_partial / (A1_partial * (1 - A1_partial))^2
  varlogitA2_partial <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_partial / (A2_partial * (1 - A2_partial))^2

  grad <- matrix(0, nrow = 2, ncol = 4)
  grad[1, 1] <- integrate(function(x) dnorm(a1 + b1 * qnorm(x)), lower = e1, upper = e2)$value
  grad[1, 2] <- integrate(function(x) dnorm(a1 + b1 * qnorm(x)) * qnorm(x), lower = e1, upper = e2)$value
  grad[2, 3] <- integrate(function(x) dnorm(a2 + b2 * qnorm(x)), lower = e1, upper = e2)$value
  grad[2, 4] <- integrate(function(x) dnorm(a2 + b2 * qnorm(x)) * qnorm(x), lower = e1, upper = e2)$value

  Delta_matrix <- grad %*% estimation$cov_matrix %*% t(grad)
  covlogit <- Delta_matrix[1, 2] / (A1_partial * A2_partial * (1 - A1_partial) * (1 - A2_partial))

  Z_partial <- (log(A1_partial / (1 - A1_partial)) - log(A2_partial / (1 - A2_partial))) / sqrt(varlogitA1_partial + varlogitA2_partial - 2 * covlogit)

  partial_area_logit_result <- ifelse(abs(Z_partial) < qnorm(1 - alpha / 2),
                                      "Using logit transformation, the partial areas are not significantly different",
                                      "Using logit transformation, the partial areas are significantly different")

  # Partial area (McClish Transformation)
  Amx <- e2 - e1
  varMcClishA1 <- 4 * Amx^2 * area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_partial / (Amx^2 - A1_partial^2)^2
  varMcClishA2 <- 4 * Amx^2 * area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_partial / (Amx^2 - A2_partial^2)^2

  covMcClish <- 4 * Amx^2 * Delta_matrix[1, 2] / ((Amx^2 - A1_partial^2) * (Amx^2 - A2_partial^2))

  Z_partial_McClish <- (log((Amx + A1_partial) / (Amx - A1_partial)) - log((Amx + A2_partial) / (Amx - A2_partial))) / sqrt(varMcClishA1 + varMcClishA2 - 2 * covMcClish)

  partial_area_mcclish_result <- ifelse(abs(Z_partial_McClish) < qnorm(1 - alpha / 2),
                                        "Using McClish transformation, the partial areas are not significantly different",
                                        "Using McClish transformation, the partial areas are significantly different")

  # Return results
  result <- list(
    test_statistic_full = Z_full,
    full_area_result = full_area_result,
    test_statistic_partial_logit = Z_partial,
    partial_area_logit_result = partial_area_logit_result,
    test_statistic_partial_mcclish = Z_partial_McClish,
    partial_area_mcclish_result = partial_area_mcclish_result
  )

  return(result)
}
