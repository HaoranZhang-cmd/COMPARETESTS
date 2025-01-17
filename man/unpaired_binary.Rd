% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.1.1.R
\name{unpaired_binary}
\alias{unpaired_binary}
\title{compare the sensitivities for unpaired binary scale data}
\usage{
unpaired_binary(data1, data2, alpha)
}
\arguments{
\item{data1}{binary test results for test 1}

\item{data2}{binary test results for test 2}

\item{alpha}{significance level}
}
\value{
whether sensitivities or specificities are significantly different
#example
#data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#unpaired_binary(data1,data2,0.05)
}
\description{
compare the sensitivities for unpaired binary scale data
}
