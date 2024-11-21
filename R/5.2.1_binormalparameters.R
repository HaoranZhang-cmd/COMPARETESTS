#' test whether binormal parameters are equal for unpaired continuous data
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param alpha significance level
#' @return whether the binormal parameters are equal
#' #example
#' #data1_no<-c(136,286,281,23,200,146,220,96,100)
#' #data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' #data2_no<-c(60,126,100,40,253,46,70,17,27)
#' #data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' #roctest_continuous_unpaired(data1_no,data1_yes,data2_no,data2_yes,0.05)
#' @export
roctest_continuous_unpaired<-function(data1_no,data2_no,data1_yes,data2_yes,alpha){
  output1<-params_estimate_continuous_unpaired(data1_no,data1_yes)
  output2<-params_estimate_continuous_unpaired(data2_no,data2_yes)
  a1<-output1$a
  b1<-output1$b
  a2<-output2$a
  b2<-output2$b
  var_a1<-output1$var_a
  var_b1<-output1$var_b
  var_a2<-output2$var_a
  var_b2<-output2$var_b
  covar_a1b1<-output1$covar_ab
  covar_a2b2<-output2$covar_ab
  a12<-a1-a2
  b12<-b1-b2
  var_a12<-var_a1+var_a2
  var_b12<-var_b1+var_b2
  covar_a12b12<-covar_a1b1+covar_a2b2
  X<-(a12^2*var_b12+b12^2*var_a12-2*a12*b12*covar_a12b12)/(var_a12*var_b12-(covar_a12b12)^2)
  print(c("If the study design is unpaired,the test statistic is",X))
  if(X>qchisq(1-alpha,2)){
    print("The binormal parameters are not equal")
  }
  else{
    print("The binormal parameters are equal")
  }
}


#' test whether binormal parameters are equal for unpaired ordinal data
#' @param data1_no test results of undiseased individuals in test 1
#' @param data2_no test results of undiseased individuals in test 2
#' @param data1_yes test results of diseased individuals in test 1
#' @param data2_yes test results of diseased individuals in test 2
#' @param alpha significance level
#' @return list containing the test statistic and whether the binormal parameters are equal
#' example
#' data1_no<-c(38,25,15,19,4)
#' data1_yes<-c(1,2,3,14,42)
#' data2_no<-c(70,7,5,7,12)
#' data2_yes<-c(8,2,2,2,48)
#' roctest_ordinal_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0.05)
#' @export
roctest_ordinal_unpaired <- function(data1_no, data2_no, data1_yes, data2_yes, alpha) {
  var.matrix <- matrix(rep(0, length = 16), nrow = 4)

  # 在unpaired条件下, 协方差矩阵是分块对角矩阵
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

  X <- (a12^2 * var_b12 + b12^2 * var_a12 - 2 * a12 * b12 * covar_a12b12) /
    (var_a12 * var_b12 - covar_a12b12^2)

  result <- list(
    test_statistic = X,
    significance = ifelse(X > qchisq(1 - alpha, 2),
                          "The binormal parameters are not equal",
                          "The binormal parameters are equal")
  )

  return(result)
}





#' using Zou's method to determine whether binormal parameters are equal for continuous paired data
#' test 1 and test 2 should be implemented on the same subjects
#' @param data1_no test results of test 1 for subjects without disease
#' @param data2_no test results of test 2 for subjects without disease
#' @param data1_yes test results of test 1 for subjects with disease
#' @param data2_yes test results of test 2 for subjects with disease
#' @param alpha significance level
#' @return list containing the test statistic and whether the binormal parameters are equal
#' example
#' data1_no<-c(136,286,281,23,200,146,220,96,100)
#' data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' data2_no<-c(60,126,100,40,253,46,70,17,27)
#' data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' roctest_continuous_paired(data1_no, data1_yes, data2_no, data2_yes, 0.05)
#' @export

roctest_continuous_paired <- function(data1_no, data1_yes, data2_no, data2_yes, alpha) {
  estimation <- params_estimate_continuous_paired(data1_no, data1_yes, data2_no, data2_yes)

  a12 <- estimation$a1 - estimation$a2
  b12 <- estimation$b1 - estimation$b2

  vara12 <- estimation$cov_matrix[1, 1] + estimation$cov_matrix[3, 3] - 2 * estimation$cov_matrix[1, 3]
  varb12 <- estimation$cov_matrix[2, 2] + estimation$cov_matrix[4, 4] - 2 * estimation$cov_matrix[2, 4]

  covara12b12 <- estimation$cov_matrix[1, 2] + estimation$cov_matrix[3, 4] -
    estimation$cov_matrix[1, 4] - estimation$cov_matrix[2, 3]

  X <- (a12^2 * vara12 + b12^2 * varb12 - 2 * a12 * b12 * covara12b12) /
    (vara12 * varb12 - covara12b12^2)

  result <- list(
    test_statistic = X,
    significance = ifelse(X > qchisq(1 - alpha, 2),
                          "The binormal parameters are not equal",
                          "The binormal parameters are equal")
  )

  return(result)
}



#' test whether binormal parameters are equal for paired ordinal data
#' @param data0 a matrix representing paired test results of undiseased individuals
#' @param data1 a matrix representing paired test results of diseased individuals
#' @param alpha significance level
#' @return list containing the test statistic and whether the binormal parameters are equal
#' example
#' data0<-matrix(c(36,20,8,6,0,0,3,3,1,0,0,1,0,4,0,0,0,2,3,2,2,1,2,5,2),nrow=5,byrow=TRUE)
#' data1<-matrix(c(1,2,2,2,1,0,0,0,0,2,0,0,0,2,0,0,0,0,2,0,0,0,1,8,39),nrow=5,byrow=TRUE)
#' roctest_ordinal_paired(data0, data1, 0.05)
#' @export
roctest_ordinal_paired <- function(data0, data1, alpha) {
  estimation <- params_estimate_ordinal_paired(data0, data1)

  a12 <- estimation$a1 - estimation$a2
  b12 <- estimation$b1 - estimation$b2

  vara12 <- estimation$cov_matrix[1, 1] + estimation$cov_matrix[3, 3] - 2 * estimation$cov_matrix[1, 3]
  varb12 <- estimation$cov_matrix[2, 2] + estimation$cov_matrix[4, 4] - 2 * estimation$cov_matrix[2, 4]

  covara12b12 <- estimation$cov_matrix[1, 2] + estimation$cov_matrix[3, 4] -
    estimation$cov_matrix[1, 4] - estimation$cov_matrix[2, 3]

  X <- (a12^2 * vara12 + b12^2 * varb12 - 2 * a12 * b12 * covara12b12) /
    (vara12 * varb12 - covara12b12^2)

  result <- list(
    test_statistic = X,
    significance = ifelse(X > qchisq(1 - alpha, 2),
                          "The binormal parameters are not equal",
                          "The binormal parameters are equal")
  )

  return(result)
}









