#' @title Missing at ranom by MCMC
#'
#' @param m Starting column number of the Y observations
#' @param n Ending column number of the Y observations
#' @param treatment Variable/column name containing the Treatment observations
#' @param n.chains Number of MCMC chains
#' @param n.iter Number of MCMC iterations
#' @param dat Data set containing treatment column and repeated observations arrange by columns observations
#'
#' @return A data table listing the posterior mean and sigma results
#' @export
#' @author Atanu Bhattacharjee, Akash Pawar and Bhrigu Kumar Rajbongshi
#' @examples
#' ##
#' data(gh)
#' hdmarjg(m=1,n=3,treatment="Treatment",n.chains=2,n.iter=10,dat=gh)
#' ##
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R and OpenBUGS. CRC Press.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2012). Applied longitudinal analysis (Vol. 998). John Wiley & Sons.


hdmarjg<-function(m,n,treatment,n.chains,n.iter,dat){
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




  dat2<-data1(dat)
  colnames(dat2)<-colnames(dat)
  GH <- dat2
  N<-nrow(GH)
  y<-as.matrix(GH[,m:n])		#matrix of Y values
  c<-ncol(y)
  trt<-as.numeric(GH[,"Treatment"])
  R1<-diag(1,c,c)
  R2<-R1

  mdata<-list("y","trt", "R1","R2","N")


  MARmodel<-function()
  {  for (i in 1:N)
  { y[i,1:3]~dmnorm(mu[trt[i],1:3],tau[trt[i],1:3,1:3]) }
    for (k in 1:2)
    { for(j in 1:3)
    { mu[k,j]~dnorm(0, 0.0001) }
    }
    tau[1,1:3,1:3]~dwish(R1[,],3)  				# prior for precision matrix
    tau[2,1:3,1:3]~dwish(R2[,],3)
    sigma[1,1:3,1:3]<-inverse(tau[1,1:3,1:3])		# covariance matrix
    sigma[2,1:3,1:3]<-inverse(tau[2,1:3,1:3])
    diff<-mu[2,3]-mu[1,3]
  }



  jags.params<-c('sigma','mu','diff')

  jagsfit <- jags(data=mdata, parameters.to.save=jags.params,n.chains =n.chains,n.iter=n.iter
                  , model.file=MARmodel)





  print(jagsfit)

}


utils::globalVariables(c("mu","beta0","pow","sigma1","sigma2","sigma3","sigma2_1","sigma3_2_1"))
