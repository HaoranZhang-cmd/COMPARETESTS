% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.2.2_nonparametric.R
\name{SEDT}
\alias{SEDT}
\title{function to calculate sensitivity and decision threshold at a fixed FPR}
\usage{
SEDT(data0, data1, SP)
}
\arguments{
\item{data0}{continuous test results of nondiseased subjects}

\item{data1}{continuous test results of diseased subjects}

\item{SP}{the fixed specificity}
}
\value{
the decision threshold,the pdf estimated at this threshold,and the corresponding Se and its variance
#example
#data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#SEDT(data0,data1,0.9)
}
\description{
function to calculate sensitivity and decision threshold at a fixed FPR
}
