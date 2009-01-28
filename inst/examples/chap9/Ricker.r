###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 9.                                ##
## Discrete time models                      ##
## Chapter 9.5.1                             ##
## Bifurcations in the discrete              ##
## logistic model                            ##
###############################################

windows()
ricker  <- function(N,r) N*exp(r*(1-N)) 
rseq    <- seq(1.5,4,0.01) # sequence of r-values

plot(0,0,xlim=range(rseq),ylim=c(0,5),type="n",
     xlab="r",ylab="Nt",main="discrete logistic model")

for ( r in rseq)
 {
  N  <- runif(1)
  for (i in 1:200) N <- ricker(N,r)   # spinup 
  for (i in 1:200){N <- ricker(N,r)
               points(r,N,pch=".",cex=1.5)}
}

