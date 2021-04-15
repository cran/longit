#' Bayesian multivariate normal regression with unstructured covariance  matrix for high dimensional longitudinal data.
#'
#'
#' Multivariate normal regression with group covaraites and unstructured covariance matrix.
#' @param m Starting number of column from where repeated observations begin
#' @param n Ending number of columns till where the repeated observations ends
#' @param time Timepoint information on which repeadted observations were taken
#' @param group A categorical variable either 0 or 1. i.e. Gender - 1 male and 0 female
#' @param chains Number of MCMC chains to be performed
#' @param iter Number of iterations to be performed
#' @param data High dimensional longitudinal data
#'
#' @return mvncovarout
#' @import R2jags
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(repdata)
#' mvncovar2(m=4,n=7,time="Age",group="Gender",chains=4,iter=100,data=repdata)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.
mvncovar2 <- function(m,n,time,group,chains,iter,data){
  timepoints<-time
  gendergrp<-group
covar2mod <- function(){
  beta ~ dnorm(0.0, 0.001)
  beta2 ~ dnorm(0.0, 0.001)
  beta4 ~ dnorm(0.0,.0001)
  for (i in 1:N2) {
    Y[i, 1:M2] ~ dmnorm(mu2[], Omega2[,])
  }
  for(j in 1:M2) {
    mu2[j] <- beta + (beta2+beta4)* age[j]
  }
  Omega2[1 : M2, 1 : M2] ~ dwish(R[,], 4)
  Sigma2[1 : M2, 1 : M2] <- inverse(Omega2[,])
  for (i in 1:N1) {
    Z[i, 1:M1] ~ dmnorm(mu1[], Omega1[,])
  }
  for(j in 1:M1) {
    mu1[j] <- beta + beta2* age[j]
  }
  Omega1[1 : M1, 1 : M1] ~ dwish(R[,],4)
  Sigma1[1 : M2, 1 : M2] <- inverse(Omega1[,])
  delta<-beta2+beta4
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

covar2data <- list(M2 = length(Age), N2 = nrow(datam), M1 = length(Age),N1 = nrow(dataf),
                     Z = z,
                     Y = y,
                     age = Age,
                     R = matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0,
                                  0, 1, 0, 0, 0, 0, 1),4, 4))
inits <- function(){list(beta = 0, beta2 = 0,beta4 = 0)}
param = c("mu1","mu2","Omega1","Omega2","Sigma1","Sigma2")
v1.sim.covar2 <- jags(covar2data, inits, model.file = covar2mod,
                      parameters.to.save =  param ,
                      n.chains = chains, n.iter = iter, n.burnin = iter/2 )

covar2out<- v1.sim.covar2
return(covar2out)
}
utils::globalVariables(c("N2","M2","beta2","beta4","age","inverse","Omega2","N1","M1","Omega1","Gender"))
