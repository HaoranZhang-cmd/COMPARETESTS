#' function to transform the original data to normal data using Box-Cox transformation
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @return the transformed data and the estimated value of the parameter for the Box-Cox transformation
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.transformed(data0,data1)
#' @export
roc.transformed<-function(data0,data1){
  n0<-length(data0)
  n1<-length(data1)
  loglikelihood<-function(lambda){
    if(abs(lambda)>=1e-4){
      data0_transformed<-(data0^lambda-1)/lambda
      data1_transformed<-(data1^lambda-1)/lambda
    }
    else{
      data0_transformed<-log(data0)
      data1_transformed<-log(data1)
    }
    mean0<-mean(data0_transformed)
    var0<-var(data0_transformed)*(n0-1)/n0
    result0<-log(sqrt(var0))*(-n0)+sum(log(data0^(lambda-1))-(data0_transformed-mean0)^2/(2*var0))
    mean1<-mean(data1_transformed)
    var1<-var(data1_transformed)*(n1-1)/length(data1)
    result1<-log(sqrt(var1))*(-n1)+sum(log(data1^(lambda-1))-(data1_transformed-mean1)^2/(2*var1))
    return(-result0-result1)
  }
  lambda<-nlminb(start=1,loglikelihood)$par
  if(abs(lambda)>=1e-4){
    data0_transformed<-(data0^lambda-1)/lambda
    data1_transformed<-(data1^lambda-1)/lambda
  }
  else{
    data0_transformed<-log(data0)
    data1_transformed<-log(data1)
  }
  output<-list(data0_transformed=NULL,data1_transformed=NULL,lambda=NULL)
  output$data0_transformed<-data0_transformed
  output$data1_transformed<-data1_transformed
  output$lambda<-lambda
  return(output)
}

#' function to calculate sensitivity and decision threshold at a fixed FPR
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @param SP the fixed specificity
#' @return the decision threshold,the pdf estimated at this threshold,and the corresponding Se and its variance
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #SEDT(data0,data1,0.9)
#' @export
SEDT<-function(data0,data1,SP){
  n0<-length(data0)
  n1<-length(data1)
  TSP<-sort(data0)[ceiling(n0*SP)]
  Se.Greenhouse<-0
  for(i in 1:n1){
    if(data1[i]>TSP){
      Se.Greenhouse<-Se.Greenhouse+1/n1
    }
  }
  varSe.Greenhouse<-Se.Greenhouse*(1-Se.Greenhouse)/n1
  new_data0<-roc.transformed(data0,data1)$data0_transformed
  new_data1<-roc.transformed(data0,data1)$data1_transformed
  new_TSP<-sort(new_data0)[ceiling(n0*SP)]
  h0<-0.9*min(sd(new_data0),(sort(new_data0)[ceiling(n0*0.75)]-sort(new_data0)[ceiling(n0*0.25)])/1.34)/(n0)^(1/5)
  h1<-0.9*min(sd(new_data1),(sort(new_data1)[ceiling(n1*0.75)]-sort(new_data1)[ceiling(n1*0.25)])/1.34)/(n1)^(1/5)
  f0<-0
  f1<-0
  for(i in 1:n0){
    f0<-f0+dnorm((new_TSP-new_data0[i])/h0)/(n0*h0)
  }
  for(i in 1:n1){
    f1<-f1+dnorm((new_TSP-new_data1[i])/h1)/(n1*h1)
  }
  varTSP<-SP*(1-SP)/(f0^2*n0)
  varSe.Linnet<-varSe.Greenhouse+varTSP*f1^2

  result<-list(TSP=NULL,f0=NULL,f1=NULL,Se=NULL,varSe=NULL)
  result$TSP<-TSP
  result$f0<-f0
  result$f1<-f1
  result$Se<-Se.Greenhouse
  result$varSe<-varSe.Linnet
  return(result)
}

#' compare sensitivities at a particular point for unpaired design and continuous data
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param alpha significance level
#' @param SP specificity of this point
#' @return list containing the test statistic and whether the sensitivities are significantly different
#' #example
#' data1_no<-c(136,286,281,23,200,146,220,96,100)
#' data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no<-c(60,126,100,40,253,46,70,17,27)
#' data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' compareSe_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0.9, 0.05)
#' @export
compareSe_unpaired <- function(data1_no, data2_no, data1_yes, data2_yes, SP, alpha) {
  Se1 <- SEDT(data1_no, data1_yes, SP)$Se
  Se2 <- SEDT(data2_no, data2_yes, SP)$Se

  varSe1 <- SEDT(data1_no, data1_yes, SP)$varSe
  varSe2 <- SEDT(data2_no, data2_yes, SP)$varSe

  Z <- (Se1 - Se2) / sqrt(varSe1 + varSe2)

  result <- list(
    test_statistic = Z,
    significance = ifelse(abs(Z) > qnorm(1 - alpha / 2),
                          "The sensitivities are significantly different",
                          "The sensitivities are not significantly different")
  )

  return(result)
}


#' compare sensitivities at a particular point for paired design and continuous data
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param alpha significance level
#' @param SP specificity of this point
#' @return list containing the test statistic and whether the sensitivities are significantly different
#' example
#' data1_no<-c(136,286,281,23,200,146,220,96,100)
#' data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no<-c(60,126,100,40,253,46,70,17,27)
#' data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' compareSe_paired(data1_no, data2_no, data1_yes, data2_yes, 0.9, 0.05)
#' @export
compareSe_paired <- function(data1_no, data2_no, data1_yes, data2_yes, SP, alpha) {
  T1SP <- SEDT(data1_no, data1_yes, SP)$TSP
  T2SP <- SEDT(data2_no, data2_yes, SP)$TSP
  f1_T1SP <- SEDT(data1_no, data1_yes, SP)$f1
  f0_T1SP <- SEDT(data1_no, data1_yes, SP)$f0
  f1_T2SP <- SEDT(data2_no, data2_yes, SP)$f1
  f0_T2SP <- SEDT(data2_no, data2_yes, SP)$f0

  nodisease <- c(data1_no, data2_no)
  disease <- c(data1_yes, data2_yes)

  matrix.no <- matrix(c(0, 0, 0, 0), nrow = 2)
  matrix.yes <- matrix(c(0, 0, 0, 0), nrow = 2)

  matrix.no[1, 1] <- sum(nodisease > T1SP & nodisease > T2SP)
  matrix.no[1, 2] <- sum(nodisease > T2SP & nodisease <= T1SP)
  matrix.no[2, 1] <- sum(nodisease > T1SP & nodisease <= T2SP)
  matrix.no[2, 2] <- sum(nodisease <= T1SP & nodisease <= T2SP)

  matrix.yes[1, 1] <- sum(disease > T1SP & disease > T2SP)
  matrix.yes[1, 2] <- sum(disease > T2SP & disease <= T1SP)
  matrix.yes[2, 1] <- sum(disease > T1SP & disease <= T2SP)
  matrix.yes[2, 2] <- sum(disease <= T1SP & disease <= T2SP)

  Se1 <- SEDT(data1_no, data1_yes, SP)$Se
  Se2 <- SEDT(data2_no, data2_yes, SP)$Se

  n0 <- length(nodisease)
  n1 <- length(disease)

  var <- SP * (1 - SP) / n0 * ((f1_T1SP / f0_T1SP)^2 + (f1_T2SP / f0_T2SP)^2) -
    2 * f1_T1SP * f1_T2SP / (n0 * f0_T1SP * f0_T2SP) * (matrix.no[1, 1] / n0 - (1 - SP)^2)

  Z <- (Se1 - Se2) / sqrt(var)

  result <- list(
    test_statistic = Z,
    significance = ifelse(abs(Z) > qnorm(1 - alpha / 2),
                          "The sensitivities are significantly different",
                          "The sensitivities are not significantly different")
  )

  return(result)
}



