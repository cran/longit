#' @title Missing not at random by MCMC
#'
#' @param m Starting column number of repeated observations
#' @param n Ending column number of the repeated observations
#' @param treatment Variable/column name containing the Treatment observations
#' @param n.chains Number of MCMC chains
#' @param n.iter Number of MCMC iterations
#' @param dat Data set containing treatment column and repeated observations
#'
#' @return Results containing a data table listing the means and sigma results
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(gh)
#' hdmnarjg(m=1,n=3,treatment="Treatment",n.chains=2,n.iter=10,dat=gh)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.
hdmnarjg<-function (m, n, treatment, n.chains,n.iter, dat)
{
  data1=function(d){
    data1=matrix(ncol=ncol(d),nrow=nrow(d))
    for(i in 1:nrow(d))
    {
      for(j in 1:ncol(d)){
        data1[i,j]=ifelse(is.na(d[i,j])==TRUE,0.0001,d[i,j])
      }
    }
    data1
  }




  data<-data1(dat)
  colnames(data)<-colnames(dat)
  y <- as.matrix(data[, m:n])
  GH <- dat
  N <- nrow(GH)
  Y <- as.matrix(GH[, m:n])

  trt <- as.numeric(GH[, treatment])
  a <- matrix(1, 2, ncol(y))
  indx <- apply(Y, 2, function(x) any(is.na(x)))
  if (indx[2] == TRUE) {
    nastat <- is.na(Y[, 2])
    N1 <- length(which(nastat == TRUE))

  }else{
    N1 <- 0
  }
  if (indx[3] == TRUE) {
    nastat <- is.na(Y[, 3])
    N2 <- length(which(nastat == TRUE))

  }else {
    N2 <- 0
  }
  dat1 <- list("y","trt", "a", "N1",
               "N2","N")

  MNARmodel <- function() {
    for (i in 1:N1) {
      y[i, 1] ~ dnorm(mu1[trt[i]], tau1[trt[i]])

    }
    for (i in (N1 + 1):N2) {
      y[i, 1] ~ dnorm(mu2[trt[i]], tau2[trt[i]])
      y[i, 2] ~ dnorm(mu2_1[i], tau2_1[trt[i]])
      mu2_1[i] <- beta0[trt[i]] + beta1[trt[i]] * y[i,
                                                    1]

    }
    for (i in (N2 + 1):N) {
      y[i, 1] ~ dnorm(mu3[trt[i]], tau3[trt[i]])
      y[i, 2] ~ dnorm(mu2_1[i], tau2_1[trt[i]])
      mu2_1[i] <- beta0[trt[i]] + beta1[trt[i]] * y[i,
                                                    1]
      y[i, 3] ~ dnorm(mu3_2_1[i], tau3_2_1[trt[i]])
      mu3_2_1[i] <- beta2[trt[i]] + beta3[trt[i]] * y[i,
                                                      1] + beta4[trt[i]] * y[i, 2]

    }

    alpha[1, 1:3] ~ ddirch(a[1, 1:3])
    alpha[2, 1:3] ~ ddirch(a[2, 1:3])
    for (k in 1:2) {
      mu1[k] ~ dnorm(0, 1e-06)
      mu2[k] ~ dnorm(0, 1e-06)
      mu3[k] ~ dnorm(0, 1e-06)
      beta0[k] ~ dnorm(0, 1e-04)
      beta1[k] ~ dnorm(0, 1e-04)
      beta2[k] ~ dnorm(0, 1e-04)
      beta3[k] ~ dnorm(0, 1e-04)
      beta4[k] ~ dnorm(0, 1e-04)
      tau1[k] <- pow(sigma1[k], -2)
      tau2[k] <- pow(sigma2[k], -2)
      tau3[k] <- pow(sigma3[k], -2)
      tau2_1[k] <- pow(sigma2_1[k], -2)
      tau3_2_1[k] <- pow(sigma3_2_1[k], -2)
      sigma1[k] ~ dunif(0, 100)
      sigma2[k] ~ dunif(0, 100)
      sigma3[k] ~ dunif(0, 50)
      sigma2_1[k] ~ dunif(0, 50)
      sigma3_2_1[k] ~ dunif(0, 50)
    }
    delta0 ~ dunif(-20, 0)
    delta1 ~ dunif(-15, 0)
    for (k in 1:2) {
      beta0ss[k] <- beta0[k] + delta0
      beta2ss[k] <- beta2[k] + delta1
      mean1[k] <- alpha[k, 1] * mu1[k] + alpha[k, 2] *
        mu2[k] + alpha[k, 3] * mu3[k]
      mean2[k] <- alpha[k, 1] * (beta0ss[k] + beta1[k] *
                                   mu1[k]) + alpha[k, 2] * (beta0[k] + beta1[k] *
                                                              mu2[k]) + alpha[k, 3] * (beta0[k] + beta1[k] *
                                                                                         mu3[k])
      mean3[k] <- alpha[k, 1] * (beta2ss[k] + beta3[k] *
                                   mu1[k] + beta4[k] * (beta0ss[k] + beta1[k] *
                                                          mu1[k])) + alpha[k, 2] * (beta2ss[k] + beta3[k] *
                                                                                      mu2[k] + beta4[k] * (beta0[k] + beta1[k] * mu2[k])) +
        alpha[k, 3] * (beta2[k] + beta3[k] * mu3[k] +
                         beta4[k] * (beta0[k] + beta1[k] * mu3[k]))
    }
    diff <- mean3[2] - mean3[1]
    prob0 <- step(diff - 1e-10)
  }
  jags.params<-c("mean1", "mean2", "mean3", "diff",
                 "prob0", "sigma1", "sigma2", "sigma3",
                 "sigma2_1", "sigma3_2_1")

  jagsfit <- jags(data=dat1, parameters.to.save=jags.params,n.chains =n.chains,n.iter=n.iter,
                  model.file=MNARmodel)
  jagsfit

}

utils::globalVariables(c("delta0","delta1","alpha","mu1","mu2","mu3","alpha","beta0","delta0","delta1","mu","mu1","mu2","mu3","pow","sigma1","sigma2","sigma2_1","sigma3","sigma3_2_1"))
