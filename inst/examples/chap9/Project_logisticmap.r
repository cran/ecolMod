###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 9.                                ##
## Discrete time models                      ##
## Chapter 9.6.1 Project                     ##
## Bifurcations of the logistic map          ##
###############################################

map  <- function(x,r) r*x*(1-x) 

rseq <- seq(2.0,4,0.005) # sequence of r-values

plot(0,0,xlim=range(rseq),ylim=c(0,1.0),type="n",xlab="r",ylab="Nt",
     main="logistic map")

for ( r in rseq)
{
x  <- runif(1)            # random initial condition, in [0,1]
for (i in 1:200) x <- map(x,r)   # spinup steps
for (i in 1:200){x <- map(x,r)  ; points(r,x,pch=".",cex=1.5)}
}

