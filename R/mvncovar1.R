#' Bayesian multivariate regression with independent covariance  matrix for high dimensional longitudinal data.
#'
#'
#' Multivariate Regression with independent covariance  matrix in longitudinal datasetup with high dimensional.
#' @param m Starting number of column from where repeated observations begin
#' @param n Ending number of columns till where the repeated observations ends
#' @param time Timepoint information on which repeadted observations were taken
#' @param group A categorical variable either 0 or 1. i.e. Gender - 1 male and 0 female
#' @param chains Number of MCMC chains to be performed
#' @param iter Number of iterations to be performed
#' @param data High dimensional longitudinal data
#'
#' @return mvncovarout lists posterior omega and sigma values.
#' @import rjags
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(repdata)
#' mvncovar1(m=4,n=7,time="Age",group="Gender",chains=10,iter=100,repdata)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.

mvncovar1<-function(m,n,time,group,chains,iter,data){
  timepoints<-time
  gendergrp<-group
mvncov<-function(){
  beta1 ~ dnorm(0.0, 0.001)
  beta2 ~ dnorm(0.0, 0.001)
  beta3 ~ dnorm(0.0,.0001)
  beta4 ~ dnorm(0.0,.0001)
  for (i in 1:N2) {
    Y[i, 1:M1] ~ dmnorm(mu2[], Omega2[,])
  }
  for(j in 1:M2) {
    mu2[j] <- beta3 + (beta2+beta4)* age[j]
  }
  Omega2[1 : M2, 1 : M2] ~ dwish(R[,], 4)
  Sigma2[1 : M2, 1 : M2] <- inverse(Omega2[,])
  for (i in 1:N1) {
    Z[i, 1:M1] ~ dmnorm(mu1[], Omega1[,])
  }
  for(j in 1:M1) {
    mu1[j] <- beta1 + beta2* age[j]
  }
  Omega1[1 : M1, 1 : M1] ~ dwish(R[,], 4)
  Sigma1[1 : M2, 1 : M2] <- inverse(Omega1[,])
  d13<-beta1-beta3
}


datam <- subset(data, gendergrp != 0)
dataf <- subset(data, gendergrp != 1)
Age <- data[,timepoints]
z = as.matrix(dataf[,m:n])
y = as.matrix(datam[,m:n])
colnames(z)<-NULL
colnames(y)<-NULL
Age <- Age[1:ncol(z)]
length(Age)
mvncovardata <- list(M2 = length(Age), N2 = nrow(datam), M1 = length(Age),N1 = nrow(dataf),
     Z = z,
     Y = y,
     age = Age,
     R = matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0,
                 0, 1, 0, 0, 0, 0, 1),4, 4))
inits <- function(){list(beta1 = 0, beta2 = 0, beta3 = 0,beta4 = 0)}
param = c("mu1","mu2","Omega1","Omega2","Sigma1","Sigma2")
v1.sim.mvncov <- jags(mvncovardata, inits, model.file =mvncov,
                    parameters.to.save = param,
                    n.chains = chains, n.iter = iter, n.burnin = iter/2 )

mvncovar1out<- v1.sim.mvncov
return(mvncovar1out)
}
utils::globalVariables(c("N2","M2","beta3","beta2","beta4","age","inverse","Omega2","N1","M1","beta1","Omega1","Gender"))
