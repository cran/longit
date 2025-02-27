% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdmnarjg.R
\name{hdmnarjg}
\alias{hdmnarjg}
\title{Missing not at random by MCMC}
\usage{
hdmnarjg(m, n, treatment, n.chains, n.iter, dat)
}
\arguments{
\item{m}{Starting column number of repeated observations}

\item{n}{Ending column number of the repeated observations}

\item{treatment}{Variable/column name containing the Treatment observations}

\item{n.chains}{Number of MCMC chains}

\item{n.iter}{Number of MCMC iterations}

\item{dat}{Data set containing treatment column and repeated observations}
}
\value{
Results containing a data table listing the means and sigma results
}
\description{
Missing not at random by MCMC
}
\examples{
##
data(gh)
hdmnarjg(m=1,n=3,treatment="Treatment",n.chains=2,n.iter=10,dat=gh)
##
}
\references{
Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.

Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.
}
\author{
Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
}
