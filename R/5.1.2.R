#' compare sensitivities for clustered binary data
#' @param data1 binary test results for test 1. It contains two rows.
#'              The first row is the number of true positive elements in each cluster
#'              The second row is the number of positive elements on gold standard in each cluster
#' @param data2 binary test results for test 2. It contains two rows.
#'              The first row is the number of true positive elements in each cluster
#'              The second row is the number of positive elements on gold standard in each cluster
#' @param alpha significance level
#' @return list containing the test statistic and whether sensitivities are significantly different
#' example
#' data1<-matrix(c(0,2,2,0,2,2,1,1,0,1,2,0,1,2,1,1,1,1,1,0,1,2,1,2,0,1,2,2,1,2,2,1,1,1,1,2,1,3,2,1,1,1,2,2,2,1,2,2,2,1),nrow=2,byrow=TRUE)
#' data2<-matrix(c(1,2,2,1,2,2,1,1,1,1,2,0,2,2,1,1,1,2,1,0,1,2,2,2,0,1,2,2,1,2,2,1,1,1,1,2,1,3,2,1,1,1,2,2,2,1,2,2,2,1),nrow=2,byrow=TRUE)
#' clustered_binary_Se(data1, data2, 0.05)
#' @export
clustered_binary_Se <- function(data1, data2, alpha) {
  if(dim(data1)[1] != 2 || dim(data2)[1] != 2) {
    stop("The row for table data should be 2")
  }

  Se1 <- data1[1,] / data1[2,]
  Se2 <- data2[1,] / data2[2,]

  Se1.total <- sum(Se1 * data1[2,] / sum(data1[2,]))
  Se2.total <- sum(Se2 * data2[2,] / sum(data2[2,]))

  Se1.var <- sum((data1[2,] / mean(data1[2,]))^2 * (Se1 - Se1.total)^2) / (length(Se1) * (length(Se1) - 1))
  Se2.var <- sum((data2[2,] / mean(data2[2,]))^2 * (Se2 - Se2.total)^2) / (length(Se2) * (length(Se2) - 1))

  Cov <- sum((data1[2,] / mean(data1[2,]))^2 * (Se1 - (Se1.total + Se2.total) / 2) *
               (Se2 - (Se1.total + Se2.total) / 2)) / (length(Se1) * (length(Se1) - 1))

  Z <- (Se1.total - Se2.total) / sqrt(Se1.var + Se2.var - 2 * Cov)

  result <- list(
    test_statistic = Z,
    significance = ifelse(Z > qnorm(alpha / 2) && Z < qnorm(1 - alpha / 2),
                          "Sensitivities are not significantly different",
                          "Sensitivities are significantly different")
  )

  return(result)
}
#' compare specificities for clustered binary data
#' @param data1 binary test results for test 1. It contains two rows.
#'              The first row is the number of true negative elements in each cluster
#'              The second row is the number of negative elements on gold standard in each cluster
#' @param data2 binary test results for test 2. It contains two rows.
#'              The first row is the number of true negative elements in each cluster
#'              The second row is the number of negative elements on gold standard in each cluster
#' @param alpha significance level
#' @return list containing the test statistic and whether specificities are significantly different
#' example
#' data1<-matrix(c(0,2,2,0,2,2,1,1,0,1,2,0,1,2,1,1,1,1,1,0,1,2,1,2,0,1,2,2,1,2,2,1,1,1,1,2,1,3,2,1,1,1,2,2,2,1,2,2,2,1),nrow=2,byrow=TRUE)
#' data2<-matrix(c(1,2,2,1,2,2,1,1,1,1,2,0,2,2,1,1,1,2,1,0,1,2,2,2,0,1,2,2,1,2,2,1,1,1,1,2,1,3,2,1,1,1,2,2,2,1,2,2,2,1),nrow=2,byrow=TRUE)
#' clustered_binary_Sp(data1, data2, 0.05)
#' @export
clustered_binary_Sp <- function(data1, data2, alpha) {
  if(dim(data1)[1] != 2 || dim(data2)[1] != 2) {
    stop("The row for table data should be 2")
  }

  Sp1 <- data1[1,] / data1[2,]
  Sp2 <- data2[1,] / data2[2,]

  Sp1.total <- sum(Sp1 * data1[2,] / sum(data1[2,]))
  Sp2.total <- sum(Sp2 * data2[2,] / sum(data2[2,]))

  Sp1.var <- sum((data1[2,] / mean(data1[2,]))^2 * (Sp1 - Sp1.total)^2) / (length(Sp1) * (length(Sp1) - 1))
  Sp2.var <- sum((data2[2,] / mean(data2[2,]))^2 * (Sp2 - Sp2.total)^2) / (length(Sp2) * (length(Sp2) - 1))

  Cov <- sum((data1[2,] / mean(data1[2,]))^2 * (Sp1 - (Sp1.total + Sp2.total) / 2) *
               (Sp2 - (Sp1.total + Sp2.total) / 2)) / (length(Sp1) * (length(Sp1) - 1))

  Z <- (Sp1.total - Sp2.total) / sqrt(Sp1.var + Sp2.var - 2 * Cov)

  result <- list(
    test_statistic = Z,
    significance = ifelse(Z > qnorm(alpha / 2) && Z < qnorm(1 - alpha / 2),
                          "Specificities are not significantly different",
                          "Specificities are significantly different")
  )

  return(result)
}

