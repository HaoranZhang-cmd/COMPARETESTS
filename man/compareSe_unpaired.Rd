% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.2.2_nonparametric.R
\name{compareSe_unpaired}
\alias{compareSe_unpaired}
\title{compare sensitivities at a particular point for unpaired design and continuous data}
\usage{
compareSe_unpaired(data1_no, data2_no, data1_yes, data2_yes, SP, alpha)
}
\arguments{
\item{data1_no}{test results of undiseased individuals in test 1}

\item{data2_no}{test results of undiseased individuals in test 2}

\item{data1_yes}{test results of diseased individuals in test 1}

\item{data2_yes}{test results of diseased individuals in test 2}

\item{SP}{specificity of this point}

\item{alpha}{significance level}
}
\value{
list containing the test statistic and whether the sensitivities are significantly different
#example
data1_no<-c(136,286,281,23,200,146,220,96,100)
data1_yes<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1360,543)
data2_no<-c(60,126,100,40,253,46,70,17,27)
data2_yes<-c(323,671,350,156,1560,543,443,509,60,230,490,80,356,120,523,303,76,353,576)
compareSe_unpaired(data1_no, data2_no, data1_yes, data2_yes, 0.9, 0.05)
}
\description{
compare sensitivities at a particular point for unpaired design and continuous data
}
