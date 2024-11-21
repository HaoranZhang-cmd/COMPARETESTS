#' He and Escobar's method to compare partial area for continuous paired data
#' @param data1_no continuous test result of test 1 for undiseased subjects
#' @param data1_yes continuous test result of test 1 for diseased subjects
#' @param data2_no continuous test result of test 2 for undiseased subjects
#' @param data2_yes continuous test result of test 2 for diseased subjects
#' @param e1 lower bound of FPR for partial area
#' @param e2 upper bound of FPR for partial area
#' @param alpha significance level
#' @return list containing the test statistic and whether the partial areas are significantly different
#' example
#' data1_no <- c(136,286,281,23,200,146,220,96,100)
#' data1_yes <- c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no <- c(60,126,100,40,253,46,70,17,27)
#' data2_yes <- c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' compare_partial_He(data1_no, data1_yes, data2_no, data2_yes, 0, 1, 0.05)
#' @export
compare_partial_He <- function(data1_no, data1_yes, data2_no, data2_yes, e1, e2, alpha) {
  n0 <- length(data1_no)
  n1 <- length(data1_yes)
  #--------------partial area estimation
  partial_He1 <- 0
  partial_He2 <- 0
  for (i in 1:n1) {
    for (j in 1:n0) {
      if (data1_no[j] >= quantile(data1_no, probs = e1) & data1_no[j] <= quantile(data1_no, probs = e2)) {
        partial_He1 <- partial_He1 + Psi(data1_yes[i], data1_no[j]) / (n0 * n1)
      }
      if (data2_no[j] >= quantile(data2_no, probs = e1) & data2_no[j] <= quantile(data2_no, probs = e2)) {
        partial_He2 <- partial_He2 + Psi(data2_yes[i], data2_no[j]) / (n0 * n1)
      }
    }
  }
  #----------------variance estimation
  np <- 0
  for (j in 1:n0) {
    if (data1_no[j] >= quantile(data1_no, probs = e1) & data1_no[j] <= quantile(data1_no, probs = e2)) {
      if (data2_no[j] >= quantile(data2_no, probs = e1) & data2_no[j] <= quantile(data2_no, probs = e2)) {
        np <- np + 1
      }
    }
  }
  tau1 <- np * partial_He1 / n0
  tau2 <- np * partial_He2 / n0
  n <- n0 + n1

  # Calculate variance for the first test
  tmp1 <- 0
  tmp2 <- 0
  for (i in 1:n0) {
    if (data1_no[i] >= quantile(data1_no, probs = e1) & data1_no[i] <= quantile(data1_no, probs = e2)) {
      V10 <- 0
      for (j in 1:n1) {
        V10 <- V10 + Psi(data1_yes[j], data1_no[i]) / n0
      }
      tmp1 <- tmp1 + (V10 - tau1) ^ 2
    }
  }
  for (j in 1:n1) {
    V01 <- 0
    for (i in 1:n0) {
      if (data1_no[i] >= quantile(data1_no, probs = e1) & data1_no[i] <= quantile(data1_no, probs = e2)) {
        V01 <- V01 + Psi(data1_yes[j], data1_no[i]) / np
      }
    }
    tmp2 <- tmp2 + (V01 - tau2) ^ 2
  }
  var_He1 <- np * tmp1 / (n * n * (np - 1)) + tmp2 / (n1 * (n1 - 1))

  # Same method for the second test
  tmp1 <- 0
  tmp2 <- 0
  for (i in 1:n0) {
    if (data2_no[i] >= quantile(data2_no, probs = e1) & data2_no[i] <= quantile(data2_no, probs = e2)) {
      V10 <- 0
      for (j in 1:n1) {
        V10 <- V10 + Psi(data2_yes[j], data2_no[i]) / n0
      }
      tmp1 <- tmp1 + (V10 - tau1) ^ 2
    }
  }
  for (j in 1:n1) {
    V01 <- 0
    for (i in 1:n0) {
      if (data2_no[i] >= quantile(data2_no, probs = e1) & data2_no[i] <= quantile(data2_no, probs = e2)) {
        V01 <- V01 + Psi(data2_yes[j], data2_no[i]) / np
      }
    }
    tmp2 <- tmp2 + (V01 - tau2) ^ 2
  }
  var_He2 <- np * tmp1 / (n * n * (np - 1)) + tmp2 / (n1 * (n1 - 1))

  # Covariance estimation
  V10T10 <- rep(0, length = n0)
  V10T20 <- rep(0, length = n0)
  V01T11 <- rep(0, length = n1)
  V01T21 <- rep(0, length = n1)
  for (i in 1:n0) {
    for (j in 1:n1) {
      if (data1_no[i] >= quantile(data1_no, probs = e1) & data1_no[i] <= quantile(data1_no, probs = e2)) {
        V10T10[i] <- V10T10[i] + Psi(data1_yes[j], data1_no[i]) / n0
      }
      if (data2_no[i] >= quantile(data2_no, probs = e1) & data2_no[i] <= quantile(data2_no, probs = e2)) {
        V10T20[i] <- V10T20[i] + Psi(data2_yes[j], data2_no[i]) / n0
      }
    }
  }
  for (j in 1:n1) {
    for (i in 1:n0) {
      if (data1_no[i] >= quantile(data1_no, probs = e1) & data1_no[i] <= quantile(data1_no, probs = e2)) {
        V01T11[j] <- V01T11[j] + Psi(data1_yes[j], data1_no[i]) / np
      }
      if (data2_no[i] >= quantile(data2_no, probs = e1) & data2_no[i] <= quantile(data2_no, probs = e2)) {
        V01T21[j] <- V01T21[j] + Psi(data2_yes[j], data2_no[i]) / np
      }
    }
  }
  S01 <- 0
  S10 <- 0
  for (j in 1:n1) {
    S01 <- S01 + (V01T11[j] - tau1) * (V01T21[j] - tau2) / (n1 - 1)
  }
  for (i in 1:n0) {
    if (data1_no[i] >= quantile(data1_no, probs = e1) & data1_no[i] <= quantile(data1_no, probs = e2)) {
      if (data2_no[i] >= quantile(data2_no, probs = e1) & data2_no[i] <= quantile(data2_no, probs = e2)) {
        S10 <- S10 + (V10T10[i] - tau1) * (V10T20[i] - tau2) / (np - 1)
      }
    }
  }
  covar <- (np / n) ^ 2 * (S10 / np + S01 / n1)

  # Calculate the test statistic
  Z <- (partial_He1 - partial_He2) / sqrt(var_He1 + var_He2 - 2 * covar)

  # Return test statistic and result
  partial_area_result <- ifelse(abs(Z) < qnorm(1 - alpha / 2),
                                "The partial areas are not significantly different",
                                "The partial areas are significantly different")

  return(list(
    test_statistic = Z,
    partial_area_result = partial_area_result
  ))
}

