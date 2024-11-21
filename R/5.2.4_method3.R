#' using Zou's method to estimate parameters for continuous data under binormal assumption
#' test 1 and test 2 should be implemented on the same subjects
#' @param data1_no test results of test 1 for subjects without disease
#' @param data2_no test results of test 2 for subjects without disease
#' @param data1_yes test results of test 1 for subjects with disease
#' @param data2_yes test results of test 2 for subjects with disease
#' @return the estimated values of the parameter for the Box-Cox transformation,the binormal parameters and the estimated covariance matrix
#' #example
#' #data1_no<-c(136,286,281,23,200,146,220,96,100)
#' #data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' #data2_no<-c(60,126,100,40,253,46,70,17,27)
#' #data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' #Zou_params_estimate(data1_no,data1_yes,data2_no,data2_yes)
#' @export
Zou_params_estimate<-function(data1_no,data1_yes,data2_no,data2_yes){
  m<-length(data1_no)
  n<-length(data1_yes)
  loglikelihood<-function(lambda){
    if(abs(lambda)<=1e-4){
      x1_star<-log(data1_no)
      x2_star<-log(data2_no)
      y1_star<-log(data1_yes)
      y2_star<-log(data2_yes)
    }
    else{
      x1_star<-(data1_no^lambda-1)/lambda
      x2_star<-(data2_no^lambda-1)/lambda
      y1_star<-(data1_yes^lambda-1)/lambda
      y2_star<-(data2_yes^lambda-1)/lambda
    }
    mu1<-mean(x1_star)
    mu2<-mean(x2_star)
    v1<-mean(y1_star)
    v2<-mean(y2_star)
    rhox<-2*sum((x1_star-mean(x1_star))*(x2_star-mean(x2_star)))/(sum((x1_star-mean(x1_star))^2)+sum((x2_star-mean(x2_star))^2))
    rhoy<-2*sum((y1_star-mean(y1_star))*(y2_star-mean(y2_star)))/(sum((y1_star-mean(y1_star))^2)+sum((y2_star-mean(y2_star))^2))
    sigma<-sqrt((sum((x1_star-mean(x1_star))^2)-2*rhox*sum((x1_star-mean(x1_star))*(x2_star-mean(x2_star)))+sum((x2_star-mean(x2_star))^2))/(2*m*(1-rhox^2)))
    tau<-sqrt((sum((y1_star-mean(y1_star))^2)-2*rhoy*sum((y1_star-mean(y1_star))*(y2_star-mean(y2_star)))+sum((y2_star-mean(y2_star))^2))/(2*n*(1-rhoy^2)))
    x1_star_trans<-(x1_star-mu1)/sigma
    x2_star_trans<-(x2_star-mu2)/sigma
    y1_star_trans<-(y1_star-v1)/tau
    y2_star_trans<-(y2_star-v2)/tau
    result<-(-m)*log(1-rhox^2)/2-2*m*log(sigma)-sum(x1_star_trans^2-2*rhox*x1_star_trans*x2_star_trans+x2_star_trans^2)/(2*(1-rhox^2))+(lambda-1)*sum(log(data1_no)+log(data2_no))-n*log(1-rhoy^2)/2-2*n*log(tau)-sum(y1_star_trans^2-2*rhoy*y1_star_trans*y2_star_trans+y2_star_trans^2)/(2*(1-rhoy^2))+(lambda-1)*(sum(log(data1_yes)+log(data2_yes)))
    return(-result)
  }
  start_param<-1
  result<-nlminb(start_param,loglikelihood)
  print(result)
  lambda<-result$par
  output<-list(a1=NULL,a2=NULL,b=NULL,cov_matrix=matrix(NA,nrow=3,ncol=3))
  #计算所有参数
  if(abs(lambda)<=1e-4){
    x1_star<-log(data1_no)
    x2_star<-log(data2_no)
    y1_star<-log(data1_yes)
    y2_star<-log(data2_yes)
  }
  else{
    x1_star<-(data1_no^lambda-1)/lambda
    x2_star<-(data2_no^lambda-1)/lambda
    y1_star<-(data1_yes^lambda-1)/lambda
    y2_star<-(data2_yes^lambda-1)/lambda
  }
  mu1<-mean(x1_star)
  mu2<-mean(x2_star)
  v1<-mean(y1_star)
  v2<-mean(y2_star)
  rhox<-2*sum((x1_star-mean(x1_star))*(x2_star-mean(x2_star)))/(sum((x1_star-mean(x1_star))^2)+sum((x2_star-mean(x2_star))^2))
  rhoy<-2*sum((y1_star-mean(y1_star))*(y2_star-mean(y2_star)))/(sum((y1_star-mean(y1_star))^2)+sum((y2_star-mean(y2_star))^2))
  sigma<-sqrt((sum((x1_star-mean(x1_star))^2)-2*rhox*sum((x1_star-mean(x1_star))*(x2_star-mean(x2_star)))+sum((x2_star-mean(x2_star))^2))/(2*m*(1-rhox^2)))
  tau<-sqrt((sum((y1_star-mean(y1_star))^2)-2*rhoy*sum((y1_star-mean(y1_star))*(y2_star-mean(y2_star)))+sum((y2_star-mean(y2_star))^2))/(2*n*(1-rhoy^2)))
  a1<-(v1-mu1)/tau
  a2<-(v2-mu2)/tau
  b<-sigma/tau
  output$cov_matrix[1,1]<-1/n+(1+rhoy^2)*a1^2/(4*n)+b^2/m
  output$cov_matrix[2,2]<-1/n+(1+rhoy^2)*a2^2/(4*n)+b^2/m
  output$cov_matrix[3,3]<-((m+1)*rhoy^2+(n+1)*rhox^2)*b^2/(4*m*n)
  output$cov_matrix[1,2]<-rhoy/n+(1+rhoy^2)*a1*a2/(4*n)+rhox*b^2/m
  output$cov_matrix[1,3]<-(1+rhoy^2)*a1*b/(4*n)
  output$cov_matrix[2,3]<-(1+rhoy^2)*a2*b/(4*n)
  output$cov_matrix[2,1]<-output$cov_matrix[1,2]
  output$cov_matrix[3,1]<-output$cov_matrix[1,3]
  output$cov_matrix[3,2]<-output$cov_matrix[2,3]
  output$a1<-a1
  output$a2<-a2
  output$b<-b
  return(output)
}

#' using Zou's method to determine whether areas are significantly different
#' test 1 and test 2 should be implemented on the same subjects
#' @param data1_no test results of test 1 for subjects without disease
#' @param data2_no test results of test 2 for subjects without disease
#' @param data1_yes test results of test 1 for subjects with disease
#' @param data2_yes test results of test 2 for subjects with disease
#' @param alpha significance level
#' @return list containing the test statistic and whether the areas are significantly different
#' example
#' data1_no <- c(136, 286, 281, 23, 200, 146, 220, 96, 100)
#' data1_yes <- c(140, 1087, 230, 183, 1256, 700, 16, 800, 253, 740, 126, 153, 283, 90, 303, 193, 76, 1360, 543)
#' data2_no <- c(60, 126, 100, 40, 253, 46, 70, 17, 27)
#' data2_yes <- c(323, 671, 350, 156, 1560, 543, 443, 509, 60, 230, 490, 80, 356, 120, 523, 303, 76, 353, 576)
#' Zou_area(data1_no, data1_yes, data2_no, data2_yes, 0.05)
#' @export
Zou_area <- function(data1_no, data1_yes, data2_no, data2_yes, alpha) {
  output <- Zou_params_estimate(data1_no, data1_yes, data2_no, data2_yes)
  a1 <- output$a1
  a2 <- output$a2
  var_a1 <- output$cov_matrix[1, 1]
  var_a2 <- output$cov_matrix[2, 2]
  covar_a1a2 <- output$cov_matrix[1, 2]

  Z <- (a1 - a2) / sqrt(var_a1 + var_a2 - 2 * covar_a1a2)

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
