#code downloaded from http://www.stat.cmu.edu/~jiashun/Research/software/NullandProp/ on Dec 15, 2015
.EstNull.func<-function (x,gamma=0.1)
{
 # x is a vector of z-values
 # gamma is a parameter, default is 0.1
 # output the estimated mean and standard deviation

 n = length(x)
 t = c(seq_len(1000))/200
 
 gan    = n^(-gamma)
 that   = 0 
 shat   = 0
 uhat   = 0
 epshat = 0

 phiplus   = rep(1,1000)
 phiminus  = rep(1,1000)
 dphiplus  = rep(1,1000)
 dphiminus = rep(1,1000)
 phi       = rep(1,1000)
 dphi      = rep(1,1000)

 for (i in seq_len(1000)) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
 }

 ind = min(c(seq(1000))[(phi - gan) <= 0])
 tt = t[ind]
 a  = phiplus[ind]
 b  = phiminus[ind]
 da = dphiplus[ind]
 db = dphiminus[ind]
 c  = phi[ind]

 that   = tt
 shat   = -(a*da + b*db)/(tt*c*c)
 shat   = sqrt(shat) 
 uhat   = -(da*b - db*a)/(c*c)
 epshat = 1 - c*exp((tt*shat)^2/2)

 return(musigma=list(mu=uhat,s=shat))
}
.epsest.func <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation

  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest=NULL

  for (j in seq_along(tt)) { 

    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in seq_len(101)) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}


