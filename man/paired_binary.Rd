% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.1.1.R
\name{paired_binary}
\alias{paired_binary}
\title{compare the sensitivities for paired binary scale data}
\usage{
paired_binary(data1, data2, alpha)
}
\arguments{
\item{data1}{test result of undiseased subjects}

\item{data2}{test result of diseased subjects}

\item{alpha}{significance level}
}
\value{
whether sensitivities or specificities are significantly different
#example
#data1<-matrix(c(12,7,11,71),nrow=2,byrow=TRUE)
#data2<-matrix(c(49,1,7,5),nrow=2,byrow=TRUE)
#paired_binary(data1,data2,0.05)
}
\description{
compare the sensitivities for paired binary scale data
}
