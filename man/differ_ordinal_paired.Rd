% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.2.3.R
\name{differ_ordinal_paired}
\alias{differ_ordinal_paired}
\title{Determining the Range of FPRs for which TPRs Differ for ordinal paired data}
\usage{
differ_ordinal_paired(data0, data1, alpha)
}
\arguments{
\item{data0}{a matrix representing paired test results for undiseased subjects}

\item{data1}{a matrix representing paired test results for diseased subjects}

\item{alpha}{significance level}
}
\value{
list containing the range of FPRs for which TPRs differ
example
data0<-matrix(c(36,20,8,6,0,0,3,3,1,0,0,1,0,4,0,0,0,2,3,2,2,1,2,5,2),nrow=5,byrow=TRUE)
data1<-matrix(c(1,2,2,2,1,0,0,0,0,2,0,0,0,2,0,0,0,0,2,0,0,0,1,8,39),nrow=5,byrow=TRUE)
differ_ordinal_paired(data0, data1, 0.05)
}
\description{
Determining the Range of FPRs for which TPRs Differ for ordinal paired data
}
