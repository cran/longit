#' Bayesian mixed effect model with MCMC
#'
#'
#' Bayesian mixed effect model with random intercepts and random slopes.
#' Fits using MCMC on longitudinal data set.
#' @param m Starting number of column from where repeated observations begin
#' @param n Ending number of columns till where the repeated observations ends
#' @param t Timepoint information on which repeadted observations were taken
#' @param group A categorical variable either 0 or 1. i.e. Gender - 1 male and 0 female
#' @param chains Number of MCMC chains to be performed
#' @param n.adapt Number of iterations to run in the JAGS adaptive phase.
#' @param data High dimensional longitudinal data
#'
#' @return Gives posterior means, standard deviation.
#' @import rjags
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(repdata)
#' Bysmixed(m=4,n=7,t="Age",group="Gender",chains=4,n.adapt=100,repdata)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.

Bysmixed<-function(m,n,t,group,chains,n.adapt,data){
  timepoints=t
  gendergrp=group
  Age <- data[,timepoints]
  Distance = as.matrix(data[,m:n])
  colnames(Distance)<-NULL
  Age <- Age[1:ncol(Distance)]

  mdata <- list(R=matrix(c(.001,0,0,.001),2,2), mu0=c(0,0)  ,
                Male=data[,gendergrp],Age = Age,
                Distance = Distance,N=nrow(data),mn=ncol(Distance))


  mixedmodinslp<-"
  model {
    for( i in 1:N){
      for(j in 1:mn){
        Distance[ i , j ] ~ dnorm(m[ i , j ], tau)
        m[ i , j  ]  <- (beta[1] + b[ i , 1 ]) + (beta[2] + b[ i , 2 ])*Age[ j ]  + beta[3]*Male[ i ]
      }
    }
    for( i in 1:N){
      b[ i , 1:2 ] ~ dmnorm(mu0[1:2], prec[1:2,1:2])
    }
    prec[1:2,1:2] ~ dwish(R[ 1:2 , 1:2 ],  2)
    for(i in 1:3){   beta[ i ] ~ dnorm(0,0.0001)  }
    tau ~dgamma(0.0001,0.0001)
    sigma <- sqrt(1/(tau*tau))
  }
  "
  m.spec<-textConnection(mixedmodinslp)

  #Inits
  inits<-function(){list(beta =c(0,0,0), tau=1, prec=matrix(c(.001,0,0,.001),2,2))}

  param = c("tau","sigma")
  jags.reg <- jags.model(m.spec,
                         data =  mdata,
                         n.chains=chains,
                         n.adapt=n.adapt)
  update(jags.reg, n.adapt)
  mixedmodout<-jags.samples(jags.reg,
                            param,
                            n.adapt)


  return(mixedmodout)
}
utils::globalVariables(c("b","tau","Male","group"))
