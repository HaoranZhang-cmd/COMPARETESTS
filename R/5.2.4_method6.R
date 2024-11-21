# compare full ROC area for ordinal clustered data
#' Input two matrices, each corresponding to a test.Each matrix consists of two rows per group.One group corresponds to one cluster.
#' The first row of each group corresponds to the detection results of diseased units in that cluster
#' The second row of each group corresponds to the detection results of undiseased units in that cluster
#' All are categorized into classes 1, 2,...,K
#' Note that two tests should be implemented on same subjects
#' @param data1 test result of clustered data for test 1
#' @param data2 test result of clustered data for test 2
#' @param alpha significance level
#' @return a list containing the test statistic and a conclusion on whether the areas are significantly different
#' example
#' data1 <- matrix(c(38,25,15,19,4,1,2,3,14,42,70,7,5,7,12,8,2,2,2,48), nrow=4, byrow=TRUE)
#' data2 <- matrix(c(38,25,19,15,4,3,2,1,14,42,70,5,7,7,12,2,2,2,8,48), nrow=4, byrow=TRUE)
#' cluster_area_compare(data1, data2, 0.05)
#' @export
cluster_area_compare <- function(data1, data2, alpha) {
  if (((dim(data1)[1] %% 2) != 0) || ((dim(data2)[1] %% 2) != 0)) {
    stop("Each cluster should contain 2 rows")
  }
  if (dim(data1)[1] != dim(data2)[1]) {
    stop("The rows of data from two tests should be the same")
  }

  K <- dim(data1)[2]
  I <- dim(data1)[1] / 2
  I10 <- 0
  I01 <- 0
  for (i in 1:dim(data1)[1]) {
    if (sum(data1[i, ]) != 0) {
      if (i %% 2 == 1) {
        I10 <- I10 + 1
      } else {
        I01 <- I01 + 1
      }
    }
  }

  n0.cluster <- rep(0, length = K)
  n1.cluster <- rep(0, length = K)
  for (i in 1:I) {
    for (j in 1:K) {
      n0.cluster[i] <- n0.cluster[i] + data1[2 * i - 1, j]
      n1.cluster[i] <- n1.cluster[i] + data1[2 * i, j]
    }
  }
  n0 <- sum(n0.cluster)
  n1 <- sum(n1.cluster)

  # Calculating the ROC area and variance for test 1
  Ac.cluster_test1 <- matrix(0, nrow = I, ncol = I)
  for (a in 1:I) {
    for (b in 1:I) {
      for (t in 1:K) {
        for (m in 1:K) {
          if (m < t) {
            Ac.cluster_test1[a, b] <- Ac.cluster_test1[a, b] + data1[2 * a - 1, t] * data1[2 * b, m]
          }
          if (m == t) {
            Ac.cluster_test1[a, b] <- Ac.cluster_test1[a, b] + data1[2 * a - 1, t] * 0.5 * data1[2 * b, m]
          }
        }
      }
    }
  }
  Ac_test1 <- sum(Ac.cluster_test1[, ]) / (n0 * n1)

  V10_test1 <- rep(0, length = I)
  V01_test1 <- rep(0, length = I)
  for (i in 1:I) {
    for (j in 1:K) {
      for (t in 1:I) {
        for (m in 1:K) {
          if (m < j) {
            V10_test1[i] <- V10_test1[i] + data1[2 * i, j] * data1[2 * t, m] / n0
          }
          if (m == j) {
            V10_test1[i] <- V10_test1[i] + data1[2 * i - 1, j] * 0.5 * data1[2 * t, m] / n0
          }
        }
      }
    }
  }

  for (i in 1:I) {
    for (j in 1:K) {
      for (t in 1:I) {
        for (m in 1:K) {
          if (m > j) {
            V01_test1[i] <- V01_test1[i] + data1[2 * i - 1, j] * data1[2 * t - 1, m] / n1
          }
          if (m == j) {
            V01_test1[i] <- V01_test1[i] + data1[2 * i - 1, j] * 0.5 * data1[2 * t - 1, m] / n1
          }
        }
      }
    }
  }

  # Variance calculations for test 1
  S10_test1 <- 0
  S01_test1 <- 0
  S11_test1 <- 0
  for (i in 1:I) {
    if (V10_test1[i] != 0) {
      S10_test1 <- S10_test1 + I10 * (V10_test1[i] - n1.cluster[i] * Ac_test1) ^ 2 / (n1 * (I10 - 1))
    }
    if (V01_test1[i] != 0) {
      S01_test1 <- S01_test1 + I01 * (V01_test1[i] - n0.cluster[i] * Ac_test1) ^ 2 / (n0 * (I01 - 1))
    }
    S11_test1 <- S11_test1 + I * (V10_test1[i] - n1.cluster[i] * Ac_test1) * (V01_test1[i] - n0.cluster[i] * Ac_test1) / (I - 1)
  }
  Ac.var_test1 <- S10_test1 / n1 + S01_test1 / n0 + 2 * S11_test1 / (n0 * n1)

  # Repeat the above for test 2
  Ac.cluster_test2 <- matrix(0, nrow = I, ncol = I)
  for (a in 1:I) {
    for (b in 1:I) {
      for (t in 1:K) {
        for (m in 1:K) {
          if (m < t) {
            Ac.cluster_test2[a, b] <- Ac.cluster_test2[a, b] + data2[2 * a - 1, t] * data2[2 * b, m]
          }
          if (m == t) {
            Ac.cluster_test2[a, b] <- Ac.cluster_test2[a, b] + data2[2 * a - 1, t] * 0.5 * data2[2 * b, m]
          }
        }
      }
    }
  }
  Ac_test2 <- sum(Ac.cluster_test2[, ]) / (n0 * n1)

  # Same variance and covariance calculations for test 2
  V10_test2 <- rep(0, length = I)
  V01_test2 <- rep(0, length = I)
  for (i in 1:I) {
    for (j in 1:K) {
      for (t in 1:I) {
        for (m in 1:K) {
          if (m < j) {
            V10_test2[i] <- V10_test2[i] + data2[2 * i, j] * data2[2 * t, m] / n0
          }
          if (m == j) {
            V10_test2[i] <- V10_test2[i] + data2[2 * i - 1, j] * 0.5 * data2[2 * t, m] / n0
          }
        }
      }
    }
  }

  for (i in 1:I) {
    for (j in 1:K) {
      for (t in 1:I) {
        for (m in 1:K) {
          if (m > j) {
            V01_test2[i] <- V01_test2[i] + data2[2 * i - 1, j] * data2[2 * t - 1, m] / n1
          }
          if (m == j) {
            V01_test2[i] <- V01_test2[i] + data2[2 * i - 1, j] * 0.5 * data2[2 * t - 1, m] / n1
          }
        }
      }
    }
  }

  # Variance calculations for test 2
  S10_test2 <- 0
  S01_test2 <- 0
  S11_test2 <- 0
  for (i in 1:I) {
    if (V10_test2[i] != 0) {
      S10_test2 <- S10_test2 + I10 * (V10_test2[i] - n1.cluster[i] * Ac_test2) ^ 2 / (n1 * (I10 - 1))
    }
    if (V01_test2[i] != 0) {
      S01_test2 <- S01_test2 + I01 * (V01_test2[i] - n0.cluster[i] * Ac_test2) ^ 2 / (n0 * (I01 - 1))
    }
    S11_test2 <- S11_test2 + I * (V10_test2[i] - n1.cluster[i] * Ac_test2) * (V01_test2[i] - n0.cluster[i] * Ac_test2) / (I - 1)
  }
  Ac.var_test2 <- S10_test2 / n1 + S01_test2 / n0 + 2 * S11_test2 / (n0 * n1)

  # Covariance calculation
  component <- c(0, 0, 0, 0)
  for (i in 1:I) {
    if (V10_test1[i] != 0) {
      component[1] <- component[1] + I10 * (V10_test1[i] - (n1.cluster[i] * Ac_test1)) * (V10_test2[i] - (n1.cluster[i] * Ac_test2)) / ((I10 - 1) * n1)
    }
    if (V01_test1[i] != 0) {
      component[2] <- component[2] + I01 * (V01_test1[i] - (n0.cluster[i] * Ac_test1)) * (V01_test2[i] - (n0.cluster[i] * Ac_test2)) / ((I01 - 1) * n0)
    }
    component[3] <- component[3] + I * (V10_test1[i] - (n1.cluster[i] * Ac_test1)) * (V01_test2[i] - (n0.cluster[i] * Ac_test2)) / (I - 1)
    component[4] <- component[4] + I * (V10_test2[i] - (n1.cluster[i] * Ac_test2)) * (V01_test1[i] - (n0.cluster[i] * Ac_test1)) / (I - 1)
  }
  covar <- component[1] / n1 + component[2] / n0 + (component[3] + component[4]) / (n1 * n0)

  # Test statistic
  Z <- (Ac_test1 - Ac_test2) / sqrt(Ac.var_test1 + Ac.var_test2 - 2 * covar)

  # Output results
  conclusion <- if (abs(Z) < qnorm(1 - alpha / 2)) {
    "The areas are not significantly different"
  } else {
    "The areas are significantly different"
  }

  return(list(test_statistic = Z, conclusion = conclusion))
}



