###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 9.                                ##
## Discrete time models                      ##
## Chapter 9.5.1 , 9.5.2                     ##
## Bifurcations in the host-parasitoid model ##
## Attractors in the host-parasitoid model   ##
###############################################

require(shape)

windows()   ## A4
par (mfrow=c(2,2))

rH <- 2.82   # rate of increase
tS <- 100    # searching time
tH <- 1      # handling time
A  <- tS/tH  # attack rate
ks <- 30     # 1/tH*a

Parasite <- function(P_H,ks)
{
 P<-P_H[1] ;  H <- P_H[2]
 f  <- A*P/(ks+H)
 return(c(H*(1-exp(-f)),       
          H * exp(rH*(1-H)-f)))
}

out <- matrix(nrow=50,ncol=2)

plottraject<-function(ks)
{
P_H <- c(0.5,0.5)
for (i in 1:100) P_H <-Parasite(P_H,ks)
for (i in 1:50) {P_H <-Parasite(P_H,ks); out[i,]<-P_H }

plot (out[,1],type="l",ylim=range(out),lwd=2,xlab="t",ylab="Population",
      main=paste("ks=",ks))
lines(out[,2],lty=2)
}

#plottraject(35)


plottraject(25)
writelabel("A")
plottraject(20)
writelabel("B")
legend("topright",c("Parasitoid","Host"),lty=c(1,2),lwd=c(2,1))

ksSeq <- seq(15,35,0.05) # sequence of a-values
plot(0,0,xlim=range(ksSeq),ylim=c(0.,2),xlab="ks",ylab="Nt",main="Bifurcation diagram")

for ( ks in ksSeq)
{
P_H <- c(0.5,0.5)
for (i in 1:500)  P_H <-Parasite(P_H,ks)   # spinup steps
for (i in 1:300)  {P_H <-Parasite(P_H,ks); points(ks,P_H[2],pch=".",cex=1.5)}
}

writelabel("C")

# domain of attraction
#a    <- 0.0433
ks   <- 23.09 
dz   <- 0.0025 
xlim <- c(0.001,0.5)
ylim <- c(0.001,0.5)

Initial <- expand.grid(P = seq(xlim[1],xlim[2],dz),
                       H = seq(ylim[1],ylim[2],dz))
plot(0,0,xlim=xlim,ylim=ylim,ylab="Parasitoid initial",xlab="Host initial",
     type="n",main="Domain of attraction")

PP   <- vector(length=100)

for ( ii in 1:nrow(Initial))
{
ini <- Initial[ii,]
P_H <- unlist(ini)
for (i in 1:100) P_H<-Parasite (P_H,ks)
for (i in 1:20) {P_H <-Parasite(P_H,ks); PP[i] <- P_H[1]}

Freq <- length(unique(trunc(PP*10)))
ifelse (Freq == 4,col<-"black",col<-"white")
rect(ini$P-dz/2,ini$H-dz/2,ini$P+dz/2,ini$H+dz/2,col=col,border=col)
}

writelabel("D")
             