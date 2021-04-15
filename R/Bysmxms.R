#' Bayesian mixed model with random intercepts and random slopes for high dimensional longitudinal data
#'
#'
#' Bayesian mixed effect model with random intercepts and slopes with longitudinally measured missing data.
#' Fits using MCMC on longitudinal data set
#' @param m Starting number of column from where repeated observations begin
#' @param n Ending number of columns till where the repeated observations ends
#' @param time Timepoint information on which repeadted observations were taken
#' @param group A categorical variable either 0 or 1. i.e. Gender - 1 male and 0 female
#' @param chains Number of MCMC chains to be performed
#' @param n.adapt Number of iterations to run in the JAGS adaptive phase.
#' @param data High dimensional longitudinal data
#'
#' @return Gives posterior means, standard deviation.
#' @import AICcmodavg
#' @import missForest
#' @import rjags
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(mesrep)
#' Bysmxms(m=4,n=7,time="Age",group="Gender",chains=4,n.adapt=100,data=msrep)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.

Bysmxms<-function(m,n,time,group,chains,n.adapt,data){

  timepoints<-time
  gendergrp<-group
  data1<-data[,m:n]
  data2<-data[,1:(m-1)]

  data1.imp <- missForest(data1) #imputed tmc longit
  data1 <- data.frame(data1.imp$ximp)
  data <- data.frame(data2,data1)
  var1 <- data[,timepoints]

  Distance = as.matrix(data[,m:n])
  colnames(Distance)<-NULL
  mn=length(m:n)
  var1 <- var1[1:mn]
  var2 <- data[,gendergrp]
  mdata <- list(R=matrix(c(.001,0,0,.001),2,2), mu0=c(0,0)  ,
                a2=var2,a1 = var1,
                Distance = Distance,N=nrow(data),MN=mn)

  mixedmodinslp<-"
  model {
    for( i in 1:N){
      for(j in 1:MN){
        Distance[ i , j ] ~ dnorm(m[ i , j ], tau)
        m[ i , j  ]  <- (beta[1] + b[ i , 1 ]) + (beta[2] + b[ i , 2 ])*a1[ j ]  + beta[3]*a2[ i ]
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


  mixedmarout<- list(mixedmodout)
  return(mixedmarout)
}
utils::globalVariables(c("b","tau","Male"))
