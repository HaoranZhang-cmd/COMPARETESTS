#----Delong's method for full area
#' this function has the same definition as the function Psi in the book.It is not important.
#' @param m the same definition as the function Psi in the book.
#' @param n the same definition as the function Psi in the book.
#' @export
Psi<-function(m,n){
  if(m>n){
    return(1)
  }
  if(m<n){
    return(0)
  }
  return(1/2)
}

#' transform ordinal test results into a vector that contains test result of each individual.It is not important.
#' @param data ordinal test results
#' @export
generate<-function(data){
  result<-c()
  for(i in 1:length(data)){
    result<-c(result,rep(i,data[i]))
  }
  return(result)
}

#' using Delong's method to compute AUC and its variance
#' @param data0 ordinal test result of undiseased subject
#' @param data1 ordinal test result of diseased subjects
#' @return estimation of AUC and its variance
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #component(data0,data1)
#' @export
component<-function(data0,data1){
  x<-generate(data0)
  y<-generate(data1)
  n0<-length(x)
  n1<-length(y)
  a<-rep(0,times=length(x))
  b<-rep(0,times=length(y))
  for(j in 1:n0){
    for(i in 1:n1)
    {
      a[j]<-a[j]+Psi(y[i],x[j])/n1
    }
  }
  for(j in 1:n1){
    for(i in 1:n0)
    {
      b[j]<-b[j]+Psi(y[j],x[i])/n0
    }
  }
  V10<-mean(b)
  V01<-mean(a)
  ADL<-(V10+V01)/2
  S10<-sum((b-ADL)^2)/(n1-1)
  S01<-sum((a-ADL)^2)/(n0-1)
  ADL.var<-1/n1*S10+1/n0*S01
  output<-list(ADL=NULL,ADL.var=NULL)
  output$ADL<-ADL
  output$ADL.var<-ADL.var
  return(output)
}

#' using Delong's method to determine whether the areas are significantly different for ordinal paired data
#' @param data1_yes ordinal test results of test 1 for diseased subjects
#' @param data1_no ordinal test results of test 1 for undiseased subjects
#' @param data2_yes ordinal test results of test 2 for diseased subjects
#' @param data2_no ordinal test results of test 2 for undiseased subjects
#' @param alpha significance level
#' @return list containing the test statistic and whether the areas are significantly different
#' example
#' data1_no <- c(38, 25, 15, 19, 4)
#' data1_yes <- c(1, 2, 3, 14, 42)
#' data2_no <- c(70, 7, 5, 7, 12)
#' data2_yes <- c(8, 2, 2, 2, 48)
#' roc_compare_Delong(data1_yes, data1_no, data2_yes, data2_no, 0.05)
#' @export
roc_compare_Delong <- function(data1_yes, data1_no, data2_yes, data2_no, alpha) {
  # Paired data
  A1DL <- component(data1_no, data1_yes)$ADL
  A2DL <- component(data2_no, data2_yes)$ADL
  varA1DL <- component(data1_no, data1_yes)$ADL.var
  varA2DL <- component(data2_no, data2_yes)$ADL.var

  data1_yes <- generate(data1_yes)
  data1_no <- generate(data1_no)
  data2_yes <- generate(data2_yes)
  data2_no <- generate(data2_no)

  n0 <- length(data1_no)
  n1 <- length(data1_yes)

  V10_T11 <- rep(0, times = n1)
  V10_T21 <- rep(0, times = n1)
  V01_T10 <- rep(0, times = n0)
  V01_T20 <- rep(0, times = n0)

  for (i in 1:n1) {
    for (j in 1:n0) {
      V10_T11[i] <- V10_T11[i] + Psi(data1_yes[i], data1_no[j]) / n0
      V10_T21[i] <- V10_T21[i] + Psi(data2_yes[i], data2_no[j]) / n0
    }
  }

  for (i in 1:n0) {
    for (j in 1:n1) {
      V01_T10[i] <- V01_T10[i] + Psi(data1_yes[j], data1_no[i]) / n1
      V01_T20[i] <- V01_T20[i] + Psi(data2_yes[j], data2_no[i]) / n1
    }
  }

  S10 <- sum((V10_T11 - A1DL) * (V10_T21 - A2DL)) / (n1 - 1)
  S01 <- sum((V01_T10 - A1DL) * (V01_T20 - A2DL)) / (n0 - 1)

  ADL_covar <- 1 / n1 * S10 + 1 / n0 * S01
  Z <- (A1DL - A2DL) / sqrt(varA1DL + varA2DL - 2 * ADL_covar)

  area_result <- ifelse(abs(Z) < qnorm(1 - alpha / 2),
                        "The areas are not significantly different",
                        "The areas are significantly different")

  # Return test statistic and result
  result <- list(
    test_statistic = Z,
    area_result = area_result
  )

  return(result)
}
