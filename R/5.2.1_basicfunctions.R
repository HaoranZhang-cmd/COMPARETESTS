#' estimate the binormal parameters for unpaired ordinal data,as well as their variances and covariances
#' @param data0 ordinal test result from undiseased individuals
#' @param data1 ordinal test result from diseased individuals
#' @return estimation of decision thresholds,binormal parameters and their covariance matrix
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #params_estimate_ordinal_unpaired(data0,data1)
#' @export
params_estimate_ordinal_unpaired<-function(data0,data1){
  K<-length(data0)
  n0<-sum(data0)
  n1<-sum(data1)
  #----------------计算对数似然函数
  loglikelihood<-function(params) {
    th<-params[1:(K-1)]
    mu0<-0
    sigma0<-1
    a<-params[K]
    b<-params[K+1]
    if(any(diff(th)<0)){
      return(Inf)
    }
    if(b<=0){
      return(Inf)
    }
    mu1<-a*sigma0/b+mu0
    sigma1<-sigma0/b
    result<-data0[1]*log(pnorm(th[1],mean=mu0,sd=sigma0))+data0[K]*log(1-pnorm(th[K-1],mean=mu0,sd=sigma0))+data1[1]*log(pnorm(th[1],mean=mu1,sd=sigma1))+data1[K]*log(1-pnorm(th[K-1],mean=mu1,sd=sigma1))
    for (i in 1:(K-2)) {
      result<-result+data0[i+1]*log(pnorm(th[i+1],mean=mu0,sd=sigma0)-pnorm(th[i],mean=mu0,sd=sigma0))+data1[i+1]*log(pnorm(th[i+1],mean=mu1,sd=sigma1)-pnorm(th[i],mean=mu1,sd=sigma1))
    }
    return(-result)
  }
  start_params<-c(seq(0,1,length=K-1),1,1)
  result<-nlminb(start_params,loglikelihood)
  print(result)
  hessian<-numDeriv::hessian(loglikelihood, result$par)
  cov_matrix<-solve(hessian)
  estimation<-list(th=NULL,a=NULL,b=NULL,cov_matrix=NULL)
  estimation$th<-result$par[1:(K-1)]
  estimation$a<-result$par[K]
  estimation$b<-result$par[K+1]
  estimation$cov_matrix<-cov_matrix[(K):(K+1),(K):(K+1)]
  return(estimation)
}


#' estimate the binormal parameters for paired ordinal data,as well as their variances and covariances
#' @param data0 a matrix representing ordinal test result of undiseased individuals
#' @param data1 a matrix representing ordinal test result of diseased individuals
#' @return estimation of decision thresholds,binormal parameters and their covariance matrix
#' #example
#' #data0<-matrix(c(36,20,8,6,0,0,3,3,1,0,0,1,0,4,0,0,0,2,3,2,2,1,2,5,2),nrow=5,byrow=TRUE)
#' #data1<-matrix(c(1,2,2,2,1,0,0,0,0,2,0,0,0,2,0,0,0,0,2,0,0,0,1,8,39),nrow=5,byrow=TRUE)
#' #params_estimate_ordinal_paired(data0,data1)
#' @export
params_estimate_ordinal_paired<-function(data0,data1){
  K<-nrow(data0)
  distribution0 <- matrix(0,nrow=K,ncol=K)
  distribution1 <- matrix(0,nrow=K,ncol=K)
  n0<-sum(data0)
  n1<-sum(data1)
  #----------------计算对数似然函数
  loglikelihood<-function(params) {
    #两个test的decision thresholds
    th1<-c(-Inf,params[1:(K-1)],Inf)
    th2<-c(-Inf,params[K:(2*K-2)],Inf)
    #_前的1/2表示是哪一个test,_后的0/1表示是否患病
    mu1_0<-0
    sigma1_0<-1
    mu2_0<-0
    sigma2_0<-1
    rho_0<-params[2*K-1]
    rho_1<-params[2*K]
    a1<-params[2*K+1]
    b1<-params[2*K+2]
    a2<-params[2*K+3]
    b2<-params[2*K+4]
    mu1_1<-a1*sigma1_0/b1+mu1_0
    mu2_1<-a2*sigma2_0/b2+mu2_0
    sigma1_1<-sigma1_0/b1
    sigma2_1<-sigma2_0/b2
    #参数要在合理范围内
    if(any(diff(th1)<0)||any(diff(th2)<0)){
      return(Inf)
    }
    if(b1<=0||b2<=0){
      return(Inf)
    }
    if(sigma1_0<=0||sigma1_1<=0){
      return(Inf)
    }
    if(sigma2_0<=0||sigma2_1<=0){
      return(Inf)
    }
    if(abs(rho_0)>=1||abs(rho_1)>=1){
      return(Inf)
    }
    #undiseased individual,diseased individual的test result都服从二元正态分布
    #对于undiseased,diseased中的每个个体,计算paired test result(X,Y)的分布
    for(i in 1:K){
      for(j in 1:K) {
        distribution0[i,j]<-pmvnorm(lower=c(th1[i],th2[j]),upper=c(th1[i+1],th2[j+1]),mean=c(mu1_0,mu2_0),sigma=matrix(c(sigma1_0^2,rho_0*sigma1_0*sigma2_0,rho_0*sigma1_0*sigma2_0,sigma2_0^2),nrow=2,byrow=TRUE))
        distribution1[i,j]<-pmvnorm(lower=c(th1[i],th2[j]),upper=c(th1[i+1],th2[j+1]),mean=c(mu1_1,mu2_1),sigma=matrix(c(sigma1_1^2,rho_1*sigma1_1*sigma2_1,rho_1*sigma1_1*sigma2_1,sigma2_1^2),nrow=2,byrow=TRUE))
      }
    }
    #计算对数似然函数(不考虑多项分布系数)
    return(-sum(data0*log(distribution0))-sum(data1*log(distribution1)))
  }
  #优化并输出结果
  start_params<-c(seq(0,1,length=K-1),seq(0,1,length=K-1),0.5,0.5,1,1,1,1)
  result<-nlminb(start_params,loglikelihood)
  print(result)
  hessian<-numDeriv::hessian(loglikelihood,result$par)
  cov_matrix<-solve(hessian)
  estimation<-list(a1=NULL,b1=NULL,a2=NULL,b2=NULL,cov_matrix=NULL)
  estimation$a1<-result$par[2*K+1]
  estimation$b1<-result$par[2*K+2]
  estimation$a2<-result$par[2*K+3]
  estimation$b2<-result$par[2*K+4]
  estimation$cov_matrix<-cov_matrix[(2*K+1):(2*K+4),(2*K+1):(2*K+4)]
  return(estimation)
}

#' Box-Cox transformation for continuous data
#' @param data0 continuous test result of undiseased subject
#' @param data1 continuous test result of diseased subject
#' @return the estimated value of the parameter for the Box-Cox transformation,
#' the binormal parameters and their covariance matrix after transformation,and the corresponding AUC
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #BoxCox(data0,data1)
#' @export
BoxCox<-function(data0,data1){
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
  result<-nlminb(start=1,loglikelihood)
  lambda<-result$par
  if(abs(lambda)>=1e-4){
    data0_transformed<-(data0^lambda-1)/lambda
    data1_transformed<-(data1^lambda-1)/lambda
  }
  else{
    data0_transformed<-log(data0)
    data1_transformed<-log(data1)
  }
  a<-(mean(data1_transformed)-mean(data0_transformed))/sqrt(var(data1_transformed)*(n1-1)/n1)
  b<-sqrt(var(data0_transformed)*(n0-1)/n0)/sqrt(var(data1_transformed)*(n1-1)/n1)
  area.full<-pnorm(a/sqrt(1+b^2))
  output<-list(lambda=NULL,a=NULL,b=NULL,area.full=NULL)
  output$lambda<-lambda
  output$a<-a
  output$b<-b
  output$area.full<-area.full
  return(output)
}

#' estimate the binormal parameters for unpaired continuous data,as well as their variances and covariances
#' @param data0 continuous test result from undiseased individuals
#' @param data1 continuous test result from diseased individuals
#' @return the estimated value of the parameter for the Box-Cox transformation,the binormal parameters,
#' their variance and covariance,and the corresponding AUC
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #params_estimate_continuous_unpaired(data0,data1)
#' @export
params_estimate_continuous_unpaired<-function(data0,data1){
  n0<-length(data0)
  n1<-length(data1)
  lambda<-BoxCox(data0,data1)$lambda
  data0_transformed<-rep(0,length=n0)
  for(i in 1:n0){
    if(abs(lambda)>=1e-4){
      data0_transformed[i]<-(data0[i]^lambda-1)/lambda
    }
    else{
      data0_transformed[i]<-log(data0[i])
    }
  }
  data1_transformed<-rep(0,length=n1)
  for(i in 1:n1){
    if(abs(lambda)>=1e-4){
      data1_transformed[i]<-(data1[i]^lambda-1)/lambda
    }
    else{
      data1_transformed[i]<-log(data1[i])
    }
  }
  a<-(mean(data1_transformed)-mean(data0_transformed))/sqrt(var(data1_transformed)*(n1-1)/n1)
  b<-sqrt(var(data0_transformed)*(n0-1)/n0)/sqrt(var(data1_transformed)*(n1-1)/n1)
  area.full<-pnorm(a/sqrt(1+b^2))
  #样本量较大时用Delta方法估计方差
  if(n0>=20&&n1>=20){
    var_a<-(n0*(a^2+2)+2*n1*b^2)/(2*n0*n1)
    var_b<-(n1+n0)*b^2/(2*n0*n1)
    covar_ab<-a*b/(2*n1)
  }
  #样本量较小时用Jackknife刀切法估计方差
  else{
    a_seq<-rep(NA,n0+n1)
    b_seq<-rep(NA,n0+n1)
    for(i in 1:n0){
      data0_new<-data0[-i]
      a_seq[i]<-BoxCox(data0_new,data1)$a
      b_seq[i]<-BoxCox(data0_new,data1)$b
    }
    for(j in 1:n1){
      data1_new<-data1[-j]
      a_seq[j+n0]<-BoxCox(data0,data1_new)$a
      b_seq[j+n0]<-BoxCox(data0,data1_new)$b
    }
    var_a<-(n0+n1-1)*sum((a_seq-mean(a_seq))^2)/(n0+n1)
    var_b<-(n0+n1-1)*sum((b_seq-mean(b_seq))^2)/(n0+n1)
    covar_ab<-(n0+n1-1)*sum((a_seq-mean(a_seq))*(b_seq-mean(b_seq)))/(n0+n1)
  }
  output<-list(lambda=NULL,a=NULL,b=NULL,var_a=NULL,var_b=NULL,covar_ab=NULL,area.full=NULL)
  output$lambda<-lambda
  output$a<-a
  output$b<-b
  output$var_a<-var_a
  output$var_b<-var_b
  output$covar_ab<-covar_ab
  output$area.full<-area.full
  return(output)
}





#' estimate the binormal parameters for paired continuous data,as well as their variances and covariances
#' 这里默认data1_no,data1_yes,data2_no,data2_yes是原始数据,需要进行Box_Cox transformation.
#' using Zou's method to estimate parameters for continuous data under binormal assumption
#' test 1 and test 2 should be implemented on the same individuals
#' @param data1_no continuous test result from undiseased individuals using test 1
#' @param data1_yes continuous test result from diseased individuals using test 1
#' @param data2_no continuous test result from undiseased individuals using test 2
#' @param data2_yes continuous test result from diseased individuals using test 2
#' @return the estimated values of the parameters for two Box-Cox transformations,the corresponding binormal parameters
#' and the estimated covariance matrix
#' #example
#' #data1_no<-c(136,286,281,23,200,146,220,96,100)
#' #data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
#' #data2_no<-c(60,126,100,40,253,46,70,17,27)
#' #data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
#' #params_estimate_continuous_paired(data1_no,data1_yes,data2_no,data2_yes)
#' @export
params_estimate_continuous_paired<-function(data1_no,data1_yes,data2_no,data2_yes){
  m<-length(data1_no)
  n<-length(data1_yes)
  loglikelihood<-function(params){
    lambda1<-params[1]
    lambda2<-params[2]
    if(abs(lambda1)>=1e-4){
      x1_star<-(data1_no^lambda1-1)/lambda1
      y1_star<-(data1_yes^lambda1-1)/lambda1
    }
    else{
      x1_star<-log(data1_no)
      y1_star<-log(data1_yes)
    }
    if(abs(lambda2)>=1e-4){
      x2_star<-(data2_no^lambda2-1)/lambda2
      y2_star<-(data2_yes^lambda2-1)/lambda2
    }else{
      x2_star<-log(data2_no)
      y2_star<-log(data2_yes)
    }
    #使用二元正态分布的极大似然估计
    mu1<-mean(x1_star)
    mu2<-mean(x2_star)
    sigma1<-sqrt(sum((x1_star-mu1)^2)/m)
    sigma2<-sqrt(sum((x2_star-mu2)^2)/m)
    rhox<-sum((x1_star-mu1)*(x2_star-mu2))/(m*sigma1*sigma2)
    v1<-mean(y1_star)
    v2<-mean(y2_star)
    tau1<-sqrt(sum((y1_star-v1)^2)/n)
    tau2<-sqrt(sum((y2_star-v2)^2)/n)
    if(sigma1<=0||sigma2<=0){
      return(Inf)
    }
    if(tau1<=0||tau2<=0){
      return(Inf)
    }
    rhoy<-sum((y1_star-v1)*(y2_star-v2))/(n*tau1*tau2)
    x1_tmp<-(x1_star-mu1)/sigma1
    x2_tmp<-(x2_star-mu2)/sigma2
    y1_tmp<-(y1_star-v1)/tau1
    y2_tmp<-(y2_star-v2)/tau2
    resultx<-(-m/2)*log(1-rhox^2)-sum(x1_tmp^2-2*rhox*x1_tmp*x2_tmp+x2_tmp^2)/(2*(1-rhox^2))+(lambda1-1)*sum(log(data1_no))+(lambda2-1)*sum(log(data2_no))-m*log(2*pi*sigma1*sigma2)
    resulty<-(-n/2)*log(1-rhoy^2)-sum(y1_tmp^2-2*rhoy*y1_tmp*y2_tmp+y2_tmp^2)/(2*(1-rhoy^2))+(lambda1-1)*sum(log(data1_yes))+(lambda2-1)*sum(log(data2_yes))-n*log(2*pi*tau1*tau2)
    return(-resultx-resulty)
  }
  start_params<-c(1,1)
  result<-nlminb(start_params,loglikelihood)
  print(result)
  estimation<-result$par
  #Using Zou's method
  lambda1<-estimation[1]
  lambda2<-estimation[2]
  if(lambda1==0){
    x1_star<-log(data1_no)
    y1_star<-log(data1_yes)
  }
  else{
    x1_star<-(data1_no^lambda1-1)/lambda1
    y1_star<-(data1_yes^lambda1-1)/lambda1
  }
  if(lambda2==0){
    x2_star<-log(data2_no)
    y2_star<-log(data2_yes)
  }
  else{
    x2_star<-(data2_no^lambda2-1)/lambda2
    y2_star<-(data2_yes^lambda2-1)/lambda2
  }
  mu1<-mean(x1_star)
  mu2<-mean(x2_star)
  sigma1<-sqrt(sum((x1_star-mu1)^2)/m)
  sigma2<-sqrt(sum((x2_star-mu2)^2)/m)
  rhox<-sum((x1_star-mu1)*(x2_star-mu2))/(m*sigma1*sigma2)
  v1<-mean(y1_star)
  v2<-mean(y2_star)
  tau1<-sqrt(sum((y1_star-v1)^2)/n)
  tau2<-sqrt(sum((y2_star-v2)^2)/n)
  rhoy<-sum((y1_star-v1)*(y2_star-v2))/(n*tau1*tau2)
  a1<-(v1-mu1)/tau1
  a2<-(v2-mu2)/tau2
  b1<-sigma1/tau1
  b2<-sigma2/tau2
  cov_matrix<-matrix(0,nrow=4,ncol=4)
  #原论文中的参数的协方差阵
  cov_matrix[1,1]<-1/n+a1^2/(2*n)+b1^2/m
  cov_matrix[2,2]<-(m+n)*b1^2/(2*m*n)
  cov_matrix[3,3]<-1/n+a2^2/(2*n)+b2^2/m
  cov_matrix[4,4]<-(m+n)*b2^2/(2*m*n)
  cov_matrix[1,2]<-a1*b1/(2*n)
  cov_matrix[1,3]<-rhoy/n+rhoy^2*a1*a2/(2*n)+rhox*b1*b2/m
  cov_matrix[1,4]<-rhoy^2*a1*b2/(2*n)
  cov_matrix[2,3]<-rhoy^2*a2*b1/(2*n)
  cov_matrix[2,4]<-(m*rhoy^2+n*rhox^2)*b1*b2/(2*m*n)
  cov_matrix[3,4]<-a2*b2/(2*n)
  cov_matrix[2,1]<-cov_matrix[1,2]
  cov_matrix[3,1]<-cov_matrix[1,3]
  cov_matrix[3,2]<-cov_matrix[2,3]
  cov_matrix[4,1]<-cov_matrix[1,4]
  cov_matrix[4,2]<-cov_matrix[2,4]
  cov_matrix[4,3]<-cov_matrix[3,4]
  output<-list(a1=NULL,b1=NULL,a2=NULL,b2=NULL,cov_matrix=NULL)
  output$a1<-a1
  output$b1<-b1
  output$a2<-a2
  output$b2<-b2
  output$cov_matrix<-cov_matrix
  return(output)
}

