#' Bayesian multivariate regression with unstructured covariance  matrix for high dimensional longitudinal data.
#'
#'
#' Multivariate Regression with unstructured covariance  matrix in longitudinal datasetup with high dimensional.
#' @param m Starting number of column from where repeated observations begin
#' @param n Ending number of columns till where the repeated observations ends
#' @param chains Number of MCMC chains to be performed
#' @param n.adapt Number of iterations to be performed
#' @param data High dimensional longitudinal data
#'
#' @return Results of posterior means and standard deviation.
#' @import rjags
#' @importFrom stats step
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(repdata)
#' creg(m=4,n=7,chains=4,n.adapt=100,data=repdata)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.
creg <- function(m,n,chains,n.adapt,data){
  Age <- data$Age
  Distance = as.matrix(data[, m:n])
  Distance = as.matrix(data[, m:n])
  colnames(Distance) <- NULL
  Age <- Age[1:ncol(Distance)]
  R = matrix(c(0.001, 0, 0, 0, 0, 0.001, 0, 0, 0, 0, 0.001,
               0, 0, 0, 0, 0.001), 4, 4)
  inits = list(beta = c(0, 0, 0), precision = matrix(c(0.001,
                                                       0, 0, 0, 0, 0.001, 0, 0, 0, 0, 0.001, 0, 0, 0, 0, 0.001),
                                                     4, 4))
  mdata <- list(Age = Age, R = R, Distance = Distance)
  mvnreg <-"
 model{
    for(i in 1:16){
      Distance[ i, 1:4 ] ~ dmnorm(m[ i, 1:4], precision[1:4,1:4])
      for(j in 1:4){  m[ i, j]  <- (beta[1] + beta[3]) + (beta[2] + beta[4])*Age[ j ] }   ## Male
    }
    for( i in 17:27){
      Distance[ i, 1:4 ] ~ dmnorm(m[ i, 1:4], precision[1:4,1:4])
      for(j in 1:4){  m[ i, j]  <- beta[1] + beta[2]*Age[ j ]    }     ## Female
    }
    covariance[1:4,1:4] <- inverse(precision[1:4,1:4])
    precision[1:4,1:4] ~ dwish(R[1:4 , 1:4 ], 4)
    for(i in 1:4){   beta[ i ] ~ dnorm(0,0.0001)  }
    for(i in 1:4){
      meanB[i] <- beta[1] + beta[3] + (beta[2] + beta[4])*Age[i]
      meanG[i] <- beta[1] + beta[2]*Age[i]
      diff[i] <- meanB[i] - meanG[i]
      testdiff[i] <- step(diff[i])
    }
    rateB <- beta[2]+beta[4]
    rateG <- beta[2]
    test <- step(beta[4])
  }
"

  mvnreg.spec<-textConnection(mvnreg)

  jags.reg <- jags.model(mvnreg.spec,
                         data =  mdata,
                         n.chains=chains,
                         n.adapt=n.adapt)
  update(jags.reg, n.adapt)
  reg.sam<-jags.samples(jags.reg,
                        c('rateB','rateG','test'),
                        n.adapt)
  return(reg.sam)
}
utils::globalVariables(c("inverse","precision","step","jags","update"))


