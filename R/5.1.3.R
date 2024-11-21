#---------5.1.3 Predictive Probability of a Positive or Negative
#' compare PPVs using Leisenring's method
#' @param data1 data from nondiseased individuals in the form of table 5.2 of the textbook
#' @param data2 data from diseased individuals in the form of table 5.2 of the textbook
#' @param alpha significance level
#' @return list containing the test statistic, p-value, and whether PPVs are significantly different
#' example
#' data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#' data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#' PPV_Leisenring(data1, data2, 0.05)
#' @export
PPV_Leisenring <- function(data1, data2, alpha) {
  if (dim(data1)[1] != 2 || dim(data2)[1] != 2) {
    stop("The row for table data should be 2")
  }

  a <- data1[1, 1]
  b <- data1[1, 2]
  c <- data1[2, 1]
  d <- data1[2, 2]
  e <- data2[1, 1]
  f <- data2[1, 2]
  g <- data2[2, 1]
  h <- data2[2, 2]

  Z <- (a + b + e + f) / (2 * e + f + g + 2 * a + b + c)
  D <- (2 * e + f + g) / (2 * e + f + g + 2 * a + b + c)

  denominator <- D * (1 - D) * (e * (1 - 2 * Z)^2 + f * (1 - Z)^2 + g * Z^2) +
    D * (1 - D) * (a * (1 - 2 * Z)^2 + b * (1 - Z)^2 + c * Z^2)
  numerator <- (e * (1 - 2 * Z) + f * (1 - Z) - g * Z)^2

  X <- numerator / denominator
  p_value <- 1 - pchisq(X, 1)

  result <- list(
    test_statistic = X,
    p_value = p_value,
    significance = ifelse(X > qchisq(1 - alpha, 1),
                          "PPVs are significantly different",
                          "PPVs are not significantly different")
  )

  return(result)
}


#' compare NPVs using Leisenring's method
#' @param data1 data from nondiseased individuals in the form of table 5.2 of the textbook
#' @param data2 data from diseased individuals in the form of table 5.2 of the textbook
#' @param alpha significance level
#' @return list containing the test statistic, p-value, and whether NPVs are significantly different
#' example
#' data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#' data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#' NPV_Leisenring(data1, data2, 0.05)
#' @export
NPV_Leisenring <- function(data1, data2, alpha) {
  if (dim(data1)[1] != 2 || dim(data2)[1] != 2) {
    stop("The row for table data should be 2")
  }

  a <- data1[1, 1]
  b <- data1[1, 2]
  c <- data1[2, 1]
  d <- data1[2, 2]
  e <- data2[1, 1]
  f <- data2[1, 2]
  g <- data2[2, 1]
  h <- data2[2, 2]

  Z <- (c + d + g + h) / (2 * d + b + c + 2 * h + f + g)
  D <- (2 * d + b + c) / (2 * d + b + c + 2 * h + f + g)

  numerator <- (d * (1 - 2 * Z) + c * (1 - Z) - b * Z)^2
  denominator <- (1 - D)^2 * (d * (1 - 2 * Z)^2 + c * (1 - Z)^2 + b * Z^2) +
    D^2 * (h * (1 - 2 * Z)^2 + g * (1 - Z)^2 + f * Z^2)

  X <- numerator / denominator
  p_value <- 1 - pchisq(X, 1)

  result <- list(
    test_statistic = X,
    p_value = p_value,
    significance = ifelse(X > qchisq(1 - alpha, 1),
                          "NPVs are significantly different",
                          "NPVs are not significantly different")
  )

  return(result)
}


#' compare PPVs using Pepe's method
#' @param data1 data from nondiseased individuals in the form of table 5.2 of the textbook
#' @param data2 data from diseased individuals in the form of table 5.2 of the textbook
#' @param alpha significance level
#' @return list containing the test statistic, p-value, and whether PPVs are significantly different
#' example
#' data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#' data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#' PPV_Pepe(data1, data2, 0.05)
#' @export
PPV_Pepe <- function(data1, data2, alpha) {
  if (dim(data1)[1] != 2 || dim(data2)[1] != 2) {
    stop("The row for table data should be 2")
  }

  a <- data1[1, 1]
  b <- data1[1, 2]
  c <- data1[2, 1]
  d <- data1[2, 2]
  e <- data2[1, 1]
  f <- data2[1, 2]
  g <- data2[2, 1]
  h <- data2[2, 2]

  PPV1 <- (e + f) / (a + b + e + f)
  PPV2 <- (e + g) / (a + c + e + g)

  n0 <- sum(data1[, ])
  n1 <- sum(data2[, ])
  N <- n0 + n1

  var <- N^2 / (e + f) / (e + g) * (2 * (c + g) / N * PPV2^2 + f / N * (1 - PPV2) + g / N * (1 - 3 * PPV2))

  Z <- log(PPV1 / PPV2) / sqrt(var / N)
  p_value <- 2 * (1 - pnorm(abs(Z)))

  result <- list(
    PPV1 = PPV1,
    PPV2 = PPV2,
    test_statistic = Z,
    p_value = p_value,
    significance = ifelse(abs(Z) > qnorm(1 - alpha / 2),
                          "PPVs are significantly different",
                          "PPVs are not significantly different")
  )

  return(result)
}


#' compare NPVs using Pepe's method
#' @param data1 data from nondiseased individuals in the form of table 5.2 of the textbook
#' @param data2 data from diseased individuals in the form of table 5.2 of the textbook
#' @param alpha significance level
#' @return list containing the test statistic, p-value, and whether NPVs are significantly different
#' example
#' data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#' data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#' NPV_Pepe(data1, data2, 0.05)
#' @export
NPV_Pepe <- function(data1, data2, alpha) {
  if (dim(data1)[1] != 2 || dim(data2)[1] != 2) {
    stop("The row for table data should be 2")
  }

  a <- data1[1, 1]
  b <- data1[1, 2]
  c <- data1[2, 1]
  d <- data1[2, 2]
  e <- data2[1, 1]
  f <- data2[1, 2]
  g <- data2[2, 1]
  h <- data2[2, 2]

  NPV1 <- (c + d) / (c + d + g + h)
  NPV2 <- (b + d) / (b + d + f + h)

  n0 <- sum(data1[, ])
  n1 <- sum(data2[, ])
  N <- n0 + n1

  var <- N^2 / (b + d) / (c + d) * ((2 * d - b - c) / N * NPV2 + (b + c) / N - 2 * (d + h) * NPV2^2 / N)

  Z <- log(NPV1 / NPV2) / sqrt(var / N)
  p_value <- 2 * (1 - pnorm(abs(Z)))

  result <- list(
    NPV1 = NPV1,
    NPV2 = NPV2,
    test_statistic = Z,
    p_value = p_value,
    significance = ifelse(abs(Z) > qnorm(1 - alpha / 2),
                          "NPVs are significantly different",
                          "NPVs are not significantly different")
  )

  return(result)
}




