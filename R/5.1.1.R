#' compare the sensitivities for unpaired binary scale data
#' @param data1 binary test results for test 1
#' @param data2 binary test results for test 2
#' @param alpha significance level
#' @return whether sensitivities or specificities are significantly different
#' #example
#' #data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#' #data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#' #unpaired_binary(data1,data2,0.05)
#' @export
#-------------5.1.1 Sensitivity and Specificity
#Unpaired binary scale data
unpaired_binary<-function(data1,data2,alpha){
  if( dim(data1)[1]!=2 || dim(data1)[2]!=2 ){
    stop("The row and column for table data should be 2")
  }
  if( dim(data2)[1]!=2 || dim(data2)[2]!=2 ){
    stop("The row and column for table data should be 2")
  }
  #灵敏度
  Se1<-data1[1,1]/sum(data1[1,])
  Se2<-data2[1,1]/sum(data2[1,])
  data<-matrix(c(data1[1,1],data1[1,2],data2[1,1],data2[1,2]),nrow=2,byrow=TRUE)
  if(any(data<=40)){
    mylist<-list(p_Se=c("p-value is",fisher.test(data)$p.value),char_Se=NULL)
    if(fisher.test(data)$p.value>alpha){
      mylist$char<-c("The sensitivities are not significantly different")
    }
    else{
      mylist$char<-c("The sensitivities are significantly different")
    }
  }
  else{
    #test statistic Z
    mylist<-list(p_Se=c("p-value is",2*(1-qnorm(abs((Se1-Se2)/sqrt(Se1*(1-Se1)/sum(data1[1,])+Se2*(1-Se2)/sum(data2[1,]))))),char_Se=NULL))
    if(2*(1-qnorm(abs((Se1-Se2)/sqrt(Se1*(1-Se1)/sum(data1[1,])+Se2*(1-Se2)/sum(data2[1,])))))>=alpha){
      mylist$char<-c("The sensitivities are not significantly different")
    }
    else{
      mylist$char<-c("The sensitivities are significantly different")
    }
  }

  #特异度
  Sp1<-data1[2,2]/sum(data1[2,])
  Sp2<-data2[2,2]/sum(data2[2,])
  data<-matrix(c(data1[2,2],data1[2,1],data2[2,2],data2[2,1]),nrow=2,byrow=TRUE)
  if(any(data<=40)){
    mylist_Sp<-list(p_Sp=c("p-value is",fisher.test(data)$p.value),char_Sp=NULL)
    if(fisher.test(data)$p.value>alpha){
      mylist_Sp$char<-c("The specificities are not significantly different")
    }
    else{
      mylist_Sp$char<-c("The specificities are significantly different")
    }
  }
  else{
    mylist_Sp<-list(p_Sp=c("p-value is",2*(1-qnorm(abs((Sp1-Sp2)/sqrt(Sp1*(1-Sp1)/sum(data1[1,])+Sp2*(1-Sp2)/sum(data2[1,]))))),char_Sp=NULL))
    if(2*(1-qnorm(abs((Se1-Se2)/sqrt(Se1*(1-Se1)/sum(data1[1,])+Se2*(1-Se2)/sum(data2[1,])))))>=alpha){
      mylist_Sp$char<-c("The specificities are not significantly different")
    }
    else{
      mylist_Sp$char<-c("The specificities are significantly different")
    }
  }
  return(c(mylist,mylist_Sp))
}



#' compare the sensitivities for paired binary scale data
#' @param data1 test result of undiseased subjects
#' @param data2 test result of diseased subjects
#' @param alpha significance level
#' @return whether sensitivities or specificities are significantly different
#' #example
#' #data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#' #data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#' #paired_binary(data1,data2,0.05)
#' @export
#paired binary scale data
paired_binary<-function(data1,data2,alpha){
  mylist=list(Se_p=NULL,Se_char=NULL,Sp_p=NULL,Sp_char=NULL)
  #灵敏度
  if(dim(data1)[1]!=2 || dim(data1)[2]!=2 ){
    stop("The row and column for table data should be 2")
  }
  if( dim(data2)[1]!=2 || dim(data2)[2]!=2 ){
    stop("The row and column for table data should be 2")
  }
  if((data2[1,2]+data2[2,1])>=20){
    X<-(data2[1,2]-data2[2,1])^2/(data2[1,2]+data2[2,1])
    mylist$Se_p<-c("The test statistic for sensitivity is",X)
    if(X<=qchisq(alpha,1,lower.tail=FALSE)){
      mylist$Se_char<-c("The sensitivities are not significantly different")
    }
    else{
      mylist$Se_char<-c("The sensitivities are significantly different")
    }
  }
  else{
    p<-0
    m<-data2[1,2]+data2[2,1]
    for(j in 1:min(data2[1,2],data2[2,1])){
      p<-p+2*choose(m,j)*(1/2)^m
    }
    mylist$Se_p<-c("The estimated p-value is",p)
    if(p>alpha){
      mylist$Se_char<-c("The sensitivities are not significantly different")
    }
    else{
      mylist$Se_char<-c("The sensitivities are significantly different")
    }
  }

  #特异度
  if((data1[1,2]+data1[2,1])>=20){
    X<-(data1[1,2]-data1[2,1])^2/(data1[1,2]+data1[2,1])
    mylist$Sp_p<-c("The test statistic for specificity is",X)
    if(X<=qchisq(alpha,1,lower.tail=FALSE)){
      mylist$Sp_char<-c("The specificities are not significantly different")
    }
    else{
      mylist$Sp_char<-c("The specificities are significantly different")
    }
  }
  else{
    p<-0
    m<-data1[1,2]+data1[2,1]
    for(j in 1:min(data1[1,2],data1[2,1])){
      p<-p+2*choose(m,j)*(1/2)^m
    }
    mylist$Sp_p<-c("The estimated p-value is",p)
    if(p>alpha){
      mylist$Sp_char<-c("The specificities are not significantly different")
    }
    else{
      mylist$Sp_char<-c("The specificities are significantly different")
    }
  }
  return(mylist)
}

