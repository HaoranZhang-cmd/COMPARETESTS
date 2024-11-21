#' function to compare sensitivities at a particular point for continuous data and unpaired design
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param alpha significance level
#' @param e FPR rate of this point
#' @return list containing the test statistic and whether the sensitivities (TPRs) are significantly different
#' #example
#' data1_no<-c(136,286,281,23,200,146,220,96,100)
#' data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no<-c(60,126,100,40,253,46,70,17,27)
#' data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' TPRtest_continuous_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0.05, 0.1)
#' @export
TPRtest_continuous_unpaired <- function(data1_no, data2_no, data1_yes, data2_yes, alpha, e) {
  a1 <- params_estimate_continuous_unpaired(data1_no, data1_yes)$a
  b1 <- params_estimate_continuous_unpaired(data1_no, data1_yes)$b
  a2 <- params_estimate_continuous_unpaired(data2_no, data2_yes)$a
  b2 <- params_estimate_continuous_unpaired(data2_no, data2_yes)$b

  var_a1 <- params_estimate_continuous_unpaired(data1_no, data1_yes)$var_a
  var_b1 <- params_estimate_continuous_unpaired(data1_no, data1_yes)$var_b
  var_a2 <- params_estimate_continuous_unpaired(data2_no, data2_yes)$var_a
  var_b2 <- params_estimate_continuous_unpaired(data2_no, data2_yes)$var_b

  covar_a1b1 <- params_estimate_continuous_unpaired(data1_no, data1_yes)$covar_ab
  covar_a2b2 <- params_estimate_continuous_unpaired(data2_no, data2_yes)$covar_ab

  a12 <- a1 - a2
  b12 <- b1 - b2

  var_a12 <- var_a1 + var_a2
  var_b12 <- var_b1 + var_b2
  covar_a12b12 <- covar_a1b1 + covar_a2b2

  D <- a12 + b12 * qnorm(e)
  varD <- var_a12 + (qnorm(e))^2 * var_b12 + 2 * qnorm(e) * covar_a12b12

  Z <- D / sqrt(varD)

  result <- list(
    test_statistic = Z,
    significance = ifelse(Z > qnorm(alpha / 2) && Z < qnorm(1 - alpha / 2),
                          "The TPRs are not significantly different",
                          "The TPRs are significantly different")
  )

  return(result)
}


#' function to compare sensitivities at a particular point for continuous data and paired design
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param alpha significance level
#' @param e FPR rate of this point
#' @return list containing the test statistic and whether the sensitivities (TPRs) are significantly different
#' #example
#' data1_no<-c(136,286,281,23,200,146,220,96,100)
#' data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no<-c(60,126,100,40,253,46,70,17,27)
#' data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' TPRtest_continuous_paired(data1_no, data2_no, data1_yes, data2_yes, 0.05, 0.1)
#' @export
TPRtest_continuous_paired <- function(data1_no, data1_yes, data2_no, data2_yes, alpha, e) {
  output <- params_estimate_continuous_paired(data1_no, data1_yes, data2_no, data2_yes)

  a1 <- output$a1
  b1 <- output$b1
  a2 <- output$a2
  b2 <- output$b2

  var_a1 <- output$cov_matrix[1, 1]
  var_b1 <- output$cov_matrix[2, 2]
  var_a2 <- output$cov_matrix[3, 3]
  var_b2 <- output$cov_matrix[4, 4]

  covar_a1b1 <- output$cov_matrix[1, 2]
  covar_a2b2 <- output$cov_matrix[3, 4]
  covar_a1a2 <- output$cov_matrix[1, 3]
  covar_b1b2 <- output$cov_matrix[2, 4]
  covar_a1b2 <- output$cov_matrix[1, 4]
  covar_a2b1 <- output$cov_matrix[2, 3]

  a12 <- a1 - a2
  b12 <- b1 - b2

  var_a12 <- var_a1 + var_a2 - 2 * covar_a1a2
  var_b12 <- var_b1 + var_b2 - 2 * covar_b1b2

  covar_a12b12 <- covar_a1b1 + covar_a2b2 - covar_a1b2 - covar_a2b1

  D <- a12 + b12 * qnorm(e)
  varD <- var_a12 + (qnorm(e))^2 * var_b12 + 2 * qnorm(e) * covar_a12b12

  Z <- D / sqrt(varD)

  result <- list(
    test_statistic = Z,
    significance = ifelse(Z > qnorm(alpha / 2) && Z < qnorm(1 - alpha / 2),
                          "The TPRs are not significantly different",
                          "The TPRs are significantly different")
  )

  return(result)
}


#' function to compare sensitivities at a particular point for ordinal data and unpaired design
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param alpha significance level
#' @param e FPR rate of this point
#' @return list containing the test statistic and whether the sensitivities (TPRs) are significantly different
#' #example
#' data1_no<-c(38,25,15,19,4)
#' data1_yes<-c(1,2,3,14,42)
#' data2_no<-c(70,7,5,7,12)
#' data2_yes<-c(8,2,2,2,48)
#' TPRtest_ordinal_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0.05, 0.1)
#' @export
TPRtest_ordinal_unpaired <- function(data1_no, data2_no, data1_yes, data2_yes, alpha, e) {
  output1 <- params_estimate_ordinal_unpaired(data1_no, data1_yes)
  output2 <- params_estimate_ordinal_unpaired(data2_no, data2_yes)

  a1 <- output1$a
  a2 <- output2$a
  b1 <- output1$b
  b2 <- output2$b

  var_a1 <- output1$cov_matrix[1, 1]
  var_a2 <- output2$cov_matrix[1, 1]
  var_b1 <- output1$cov_matrix[2, 2]
  var_b2 <- output2$cov_matrix[2, 2]

  covar_a1b1 <- output1$cov_matrix[1, 2]
  covar_a2b2 <- output2$cov_matrix[1, 2]

  a12 <- a1 - a2
  b12 <- b1 - b2

  var_a12 <- var_a1 + var_a2
  var_b12 <- var_b1 + var_b2
  covar_a12b12 <- covar_a1b1 + covar_a2b2

  D <- a12 + b12 * qnorm(e)
  varD <- var_a12 + (qnorm(e))^2 * var_b12 + 2 * qnorm(e) * covar_a12b12

  Z <- D / sqrt(varD)

  result <- list(
    test_statistic = Z,
    significance = ifelse(Z > qnorm(alpha / 2) && Z < qnorm(1 - alpha / 2),
                          "The TPRs are not significantly different",
                          "The TPRs are significantly different")
  )

  return(result)
}


#' function to compare sensitivities at a particular point for ordinal data and paired design
#' @param data0 a matrix representing paired test results of undiseased individuals
#' @param data1 a matrix representing paired test results of diseased individuals
#' @param alpha significance level
#' @param e FPR rate of this point
#' @return list containing the test statistic and whether the sensitivities (TPRs) are significantly different
#' #example
#' data0<-matrix(c(36,20,8,6,0,0,3,3,1,0,0,1,0,4,0,0,0,2,3,2,2,1,2,5,2),nrow=5,byrow=TRUE)
#' data1<-matrix(c(1,2,2,2,1,0,0,0,0,2,0,0,0,2,0,0,0,0,2,0,0,0,1,8,39),nrow=5,byrow=TRUE)
#' TPRtest_ordinal_paired(data0, data1, 0.05, 0.1)
#' @export
TPRtest_ordinal_paired <- function(data0, data1, alpha, e) {
  output <- params_estimate_ordinal_paired(data0, data1)

  a1 <- output$a1
  b1 <- output$b1
  a2 <- output$a2
  b2 <- output$b2

  var_a1 <- output$cov_matrix[1, 1]
  var_b1 <- output$cov_matrix[2, 2]
  var_a2 <- output$cov_matrix[3, 3]
  var_b2 <- output$cov_matrix[4, 4]

  covar_a1b1 <- output$cov_matrix[1, 2]
  covar_a2b2 <- output$cov_matrix[3, 4]
  covar_a1a2 <- output$cov_matrix[1, 3]
  covar_b1b2 <- output$cov_matrix[2, 4]
  covar_a1b2 <- output$cov_matrix[1, 4]
  covar_a2b1 <- output$cov_matrix[2, 3]

  a12 <- a1 - a2
  b12 <- b1 - b2

  var_a12 <- var_a1 + var_a2 - 2 * covar_a1a2
  var_b12 <- var_b1 + var_b2 - 2 * covar_b1b2

  covar_a12b12 <- covar_a1b1 + covar_a2b2 - covar_a1b2 - covar_a2b1

  D <- a12 + b12 * qnorm(e)
  varD <- var_a12 + (qnorm(e))^2 * var_b12 + 2 * qnorm(e) * covar_a12b12

  Z <- D / sqrt(varD)

  result <- list(
    test_statistic = Z,
    significance = ifelse(Z > qnorm(alpha / 2) && Z < qnorm(1 - alpha / 2),
                          "The TPRs are not significantly different",
                          "The TPRs are significantly different")
  )

  return(result)
}

