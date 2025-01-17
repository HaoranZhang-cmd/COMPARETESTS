% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.2.4_method4.R
\name{component}
\alias{component}
\title{using Delong's method to compute AUC and its variance}
\usage{
component(data0, data1)
}
\arguments{
\item{data0}{ordinal test result of undiseased subject}

\item{data1}{ordinal test result of diseased subjects}
}
\value{
estimation of AUC and its variance
#example
#data0<-c(38,25,15,19,4)
#data1<-c(1,2,3,14,42)
#component(data0,data1)
}
\description{
using Delong's method to compute AUC and its variance
}
