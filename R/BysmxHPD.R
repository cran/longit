#' Bayesian mixed effect model for high dimensional longitduinal data with highest posterior density interval (HPDI).
#'
#'
#' Bayesian mixed effect model with random intercept and slopes provides inference with highest posterior density interval (HPDI).
#' Data longitudinally measured missing value and having batched
#' information.
#' Fits using MCMC on longitudinal data set
#' @param m Starting number of column from where repeated observations begin
#' @param tmax Ending number of columns till where the repeated observations ends
#' @param t Timepoint information on which repeadted observations were taken
#' @param group A categorical variable either 0 or 1. i.e. Gender - 1 male and 0 female
#' @param chains Number of MCMC chains to be performed
#' @param iter Number of iterations to be performed
#' @param out DIC/HPD outcome
#' @param data High dimensional longitudinal data
#' @return Gives posterior means, standard deviation.
#' @import R2jags
#' @import missForest
#' @import AICcmodavg
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(msrep)
#' BysmxHPD(m=c(4,8,12),tmax=4,t="Age",group="Gender",chains=4,iter=1000,out="hpD",data=msrep)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.


BysmxHPD<-function(m,tmax,t,group,chains,iter,out,data){
  #mixed model with random effects and random slopes
  timepoints=t
  gendergrp=group
  w1<-m[1]
  w2<-(m[length(m)])+tmax-1
  data1<-data[,w1:w2]
  data2<-data[,1:(w1-1)]

  data1.imp <- missForest(data1)
  data1 <- data.frame(data1.imp$ximp)
  data <- data.frame(data2,data1)
  var1 <- data[,timepoints]
  dmat <- matrix(nrow=0,ncol=2)
  rlt=list()
  for( k in m){
    m1=k
    m2=k+tmax-1
    Distance = as.matrix(data[,m1:m2])
    colnames(Distance)<-NULL
    mn=length(m1:m2)
    var1 <- var1[1:mn]
    var2 <- data[,gendergrp]
    mdata <- list(R=matrix(c(.001,0,0,.001),2,2), mu0=c(0,0)  ,
                  a2=var2,a1 = var1,
                  Distance = Distance,N=nrow(data),MN=mn)


    mixedmodinslp <- function(){
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
    inits<-function(){list(beta =c(0,0,0), tau=1, prec=matrix(c(.001,0,0,.001),2,2))}

    param = c("tau","sigma")
    v1.sim.rem1 <- jags(mdata, inits, model.file = mixedmodinslp,
                        parameters.to.save = param,
                        n.chains = chains, n.iter = iter, n.burnin = iter/2 )

    rlt[[k]] <- v1.sim.rem1

    if(out == "DIC"){
      dic.var <- DIC(v1.sim.rem1)
      var.na <- names(data[k])
      dic.data<-data.frame(var.na,dic.var)
      dmat <- rbind(dmat,dic.data)
    }


  }

  ra <- data.frame(dmat)
  if(out=="DIC"){
   colnames(ra)<-c("Variable","DIC")
    result<-ra
  }


  if(out=="hpD"){
    result<-rlt
    result = result[-which(sapply(result, is.null))]
  }
  return(result)
}
utils::globalVariables(c("b","tau","Male","MN","N","a1","a2"))
