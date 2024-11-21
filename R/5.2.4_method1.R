# This file uses equations 4.52-4.58
#' compare area and partial area for ordinal data and unpaired design
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param e1 lower bound of FPR for the partial area
#' @param e2 upper bound of FPR for the partial area
#' @param alpha significance level
#' @return list containing the test statistics and whether the areas or partial areas are significantly different
#' example
#' data1_no<-c(38,25,15,19,4)
#' data1_yes<-c(1,2,3,14,42)
#' data2_no<-c(70,7,5,7,12)
#' data2_yes<-c(8,2,2,2,48)
#' area_ordinal_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0, 0.2, 0.05)
#' @export
area_ordinal_unpaired <- function(data1_no, data2_no, data1_yes, data2_yes, e1, e2, alpha) {
  estimation1 <- params_estimate_ordinal_unpaired(data1_no, data1_yes)
  estimation2 <- params_estimate_ordinal_unpaired(data2_no, data2_yes)

  a1 <- estimation1$a
  b1 <- estimation1$b
  a2 <- estimation2$a
  b2 <- estimation2$b

  var_a1 <- estimation1$cov_matrix[1, 1]
  var_b1 <- estimation1$cov_matrix[2, 2]
  covar_a1b1 <- estimation1$cov_matrix[1, 2]

  var_a2 <- estimation2$cov_matrix[1, 1]
  var_b2 <- estimation2$cov_matrix[2, 2]
  covar_a2b2 <- estimation2$cov_matrix[1, 2]

  output1 <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)
  output2 <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)

  Z1_full <- output1$area.full
  Z2_full <- output2$area.full
  var_full <- output1$var_full + output2$var_full

  Z1_partial <- output1$area.partial
  Z2_partial <- output2$area.partial
  var_partial <- output1$var_partial + output2$var_partial

  Z_full <- (Z1_full - Z2_full) / sqrt(var_full)
  Z_partial <- (Z1_partial - Z2_partial) / sqrt(var_partial)

  result <- list(
    test_statistic_full = Z_full,
    area_significance = ifelse(abs(Z_full) < qnorm(1 - alpha / 2),
                               "The areas are not significantly different",
                               "The areas are significantly different"),
    test_statistic_partial = Z_partial,
    partial_area_significance = ifelse(abs(Z_partial) < qnorm(1 - alpha / 2),
                                       "The partial areas are not significantly different",
                                       "The partial areas are significantly different")
  )

  return(result)
}


#' compare area and partial area for continuous data and unpaired design
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param e1 lower bound of FPR for the partial area
#' @param e2 upper bound of FPR for the partial area
#' @param alpha significance level
#' @return list containing the test statistics and whether the areas or partial areas are significantly different
#' example
#' data1_no<-c(136,286,281,23,200,146,220,96,100)
#' data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no<-c(60,126,100,40,253,46,70,17,27)
#' data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' area_continuous_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0, 0.2, 0.05)
#' @export
area_continuous_unpaired <- function(data1_no, data2_no, data1_yes, data2_yes, e1, e2, alpha) {
  estimation1 <- params_estimate_continuous_unpaired(data1_no, data1_yes)
  estimation2 <- params_estimate_continuous_unpaired(data2_no, data2_yes)

  a1 <- estimation1$a
  b1 <- estimation1$b
  a2 <- estimation2$a
  b2 <- estimation2$b

  var_a1 <- estimation1$var_a
  var_b1 <- estimation1$var_b
  covar_a1b1 <- estimation1$covar_ab

  var_a2 <- estimation2$var_a
  var_b2 <- estimation2$var_b
  covar_a2b2 <- estimation2$covar_ab

  output1 <- area_estimation_continuous(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)
  output2 <- area_estimation_continuous(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)

  Z1_full <- output1$area.full
  Z2_full <- output2$area.full
  var_full <- output1$var_full + output2$var_full

  Z1_partial <- output1$area.partial
  Z2_partial <- output2$area.partial
  var_partial <- output1$var_partial + output2$var_partial

  Z_full <- (Z1_full - Z2_full) / sqrt(var_full)
  Z_partial <- (Z1_partial - Z2_partial) / sqrt(var_partial)

  result <- list(
    test_statistic_full = Z_full,
    area_significance = ifelse(abs(Z_full) < qnorm(1 - alpha / 2),
                               "The areas are not significantly different",
                               "The areas are significantly different"),
    test_statistic_partial = Z_partial,
    partial_area_significance = ifelse(abs(Z_partial) < qnorm(1 - alpha / 2),
                                       "The partial areas are not significantly different",
                                       "The partial areas are significantly different")
  )

  return(result)
}


#' compare area and partial area for continuous data and paired design
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param e1 lower bound of FPR for the partial area
#' @param e2 upper bound of FPR for the partial area
#' @param alpha significance level
#' @return list containing the test statistics and whether the areas or partial areas are significantly different
#' example
#' data1_no<-c(136,286,281,23,200,146,220,96,100)
#' data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no<-c(60,126,100,40,253,46,70,17,27)
#' data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' area_continuous_paired(data1_no, data2_no, data1_yes, data2_yes, 0, 0.2, 0.05)
#' @export
area_continuous_paired <- function(data1_no, data2_no, data1_yes, data2_yes, e1, e2, alpha) {
  output <- params_estimate_continuous_paired(data1_no, data1_yes, data2_no, data2_yes)

  a1 <- output$a1
  b1 <- output$b1
  a2 <- output$a2
  b2 <- output$b2
  a <- c(a1, a2)
  b <- c(b1, b2)
  e <- c(e1, e2)

  # 协方差
  var_a1 <- output$cov_matrix[1, 1]
  var_a2 <- output$cov_matrix[3, 3]
  var_b1 <- output$cov_matrix[2, 2]
  var_b2 <- output$cov_matrix[4, 4]
  covar_a1b1 <- output$cov_matrix[1, 2]
  covar_a2b2 <- output$cov_matrix[3, 4]
  covar_a1a2 <- output$cov_matrix[1, 3]
  covar_b1b2 <- output$cov_matrix[2, 4]
  covar_a2b1 <- output$cov_matrix[2, 3]
  covar_a1b2 <- output$cov_matrix[1, 4]

  # 其它相关参数
  h <- matrix(c(0, 0, 0, 0), nrow = 2)
  for (i in 1:2) {
    for (j in 1:2) {
      h[i, j] <- (qnorm(e[j]) + a[i] * b[i] / (1 + b[i]^2)) * sqrt(1 + b[i]^2)
    }
  }

  f_full <- c(0, 0)
  g_full <- c(0, 0)
  f_partial <- c(0, 0)
  g_partial <- c(0, 0)

  for (i in 1:2) {
    f_full[i] <- exp(-a[i]^2 / (2 * (1 + b[i]^2))) / sqrt(2 * pi * (1 + b[i]^2))
    g_full[i] <- -a[i] * b[i] * exp(-a[i]^2 / (2 * (1 + b[i]^2))) / sqrt(2 * pi * (1 + b[i]^2)^3)
    f_partial[i] <- exp(-a[i]^2 / (2 * (1 + b[i]^2))) * (pnorm(h[i, 2]) - pnorm(h[i, 1])) / sqrt(2 * pi * (1 + b[i]^2))
    g_partial[i] <- exp(-a[i]^2 / (2 * (1 + b[i]^2))) * (exp(-h[i, 1]^2 / 2) - exp(-h[i, 2]^2 / 2)) / (2 * pi * (1 + b[i]^2)) -
      a[i] * b[i] * exp(-a[i]^2 / (2 * (1 + b[i]^2))) * (pnorm(h[i, 2]) - pnorm(h[i, 1])) / sqrt(2 * pi * (1 + b[i]^2)^3)
  }

  cov_full <- f_full[1] * f_full[2] * covar_a1a2 + g_full[1] * g_full[2] * covar_b1b2 + g_full[1] * f_full[2] * covar_a2b1 + f_full[1] * g_full[2] * covar_a1b2
  cov_partial <- f_partial[1] * f_partial[2] * covar_a1a2 + g_partial[1] * g_partial[2] * covar_b1b2 + g_partial[1] * f_partial[2] * covar_a2b1 + f_partial[1] * g_partial[2] * covar_a1b2

  # 计算结果
  Z1_full <- area_estimation_continuous(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$area.full
  Z1_partial <- area_estimation_continuous(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$area.partial
  Z2_full <- area_estimation_continuous(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$area.full
  Z2_partial <- area_estimation_continuous(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$area.partial

  var_full1 <- area_estimation_continuous(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_full
  var_full2 <- area_estimation_continuous(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_full
  var_partial1 <- area_estimation_continuous(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)$var_partial
  var_partial2 <- area_estimation_continuous(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)$var_partial

  var_full <- var_full1 + var_full2 - 2 * cov_full
  var_partial <- var_partial1 + var_partial2 - 2 * cov_partial

  Z_full <- (Z1_full - Z2_full) / sqrt(var_full)
  Z_partial <- (Z1_partial - Z2_partial) / sqrt(var_partial)

  result <- list(
    test_statistic_full = Z_full,
    full_area_significance = ifelse(abs(Z_full) < qnorm(1 - alpha / 2),
                                    "The full areas are not significantly different",
                                    "The full areas are significantly different"),
    test_statistic_partial = Z_partial,
    partial_area_significance = ifelse(abs(Z_partial) < qnorm(1 - alpha / 2),
                                       "The partial areas are not significantly different",
                                       "The partial areas are significantly different")
  )

  return(result)
}


#' compare area and partial area for ordinal data and paired design
#' @param data0 a matrix representing paired test results for undiseased subjects
#' @param data1 a matrix representing paired test results for diseased subjects
#' @param e1 lower bound of FPR for the partial area
#' @param e2 upper bound of FPR for the partial area
#' @param alpha significance level
#' @return list containing the test statistics and whether the areas or partial areas are significantly different
#' example
#' data0<-matrix(c(36,20,8,6,0,0,3,3,1,0,0,1,0,4,0,0,0,2,3,2,2,1,2,5,2),nrow=5,byrow=TRUE)
#' data1<-matrix(c(1,2,2,2,1,0,0,0,0,2,0,0,0,2,0,0,0,0,2,0,0,0,1,8,39),nrow=5,byrow=TRUE)
#' area_ordinal_paired(data0, data1, 0, 0.2, 0.05)
#' @export
area_ordinal_paired <- function(data0, data1, e1, e2, alpha) {
  output <- params_estimate_ordinal_paired(data0, data1)

  a1 <- output$a1
  b1 <- output$b1
  a2 <- output$a2
  b2 <- output$b2

  # Covariances and variances
  var_a1 <- output$cov_matrix[1, 1]
  var_a2 <- output$cov_matrix[3, 3]
  var_b1 <- output$cov_matrix[2, 2]
  var_b2 <- output$cov_matrix[4, 4]

  covar_a1b1 <- output$cov_matrix[1, 2]
  covar_a2b2 <- output$cov_matrix[3, 4]
  covar_a1a2 <- output$cov_matrix[1, 3]
  covar_b1b2 <- output$cov_matrix[2, 4]
  covar_a2b1 <- output$cov_matrix[2, 3]
  covar_a1b2 <- output$cov_matrix[1, 4]

  # Estimations
  output1 <- area_estimation_ordinal(a1, b1, var_a1, var_b1, covar_a1b1, e1, e2)
  output2 <- area_estimation_ordinal(a2, b2, var_a2, var_b2, covar_a2b2, e1, e2)

  Z1_full <- output1$area.full
  Z2_full <- output2$area.full
  Z1_partial <- output1$area.partial
  Z2_partial <- output2$area.partial

  var_full1 <- output1$var_full
  var_full2 <- output2$var_full
  var_partial1 <- output1$var_partial
  var_partial2 <- output2$var_partial

  cov_full <- covar_a1a2 + covar_b1b2 + covar_a2b1 + covar_a1b2
  cov_partial <- covar_a1a2 + covar_b1b2 + covar_a2b1 + covar_a1b2

  var_full <- var_full1 + var_full2 - 2 * cov_full
  var_partial <- var_partial1 + var_partial2 - 2 * cov_partial

  Z_full <- (Z1_full - Z2_full) / sqrt(var_full)
  Z_partial <- (Z1_partial - Z2_partial) / sqrt(var_partial)

  result <- list(
    test_statistic_full = Z_full,
    full_area_significance = ifelse(abs(Z_full) < qnorm(1 - alpha / 2),
                                    "The full areas are not significantly different",
                                    "The full areas are significantly different"),
    test_statistic_partial = Z_partial,
    partial_area_significance = ifelse(abs(Z_partial) < qnorm(1 - alpha / 2),
                                       "The partial areas are not significantly different",
                                       "The partial areas are significantly different")
  )

  return(result)
}
