###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 5.                                ##
## Model solution - analytical methods       ##
## Chapter 5.5.2. Project                    ##
## Oxygen dynamics in sediment               ##
###############################################

par(mfrow=c(2,2))

###################################
# 1-dimensional oxygen model
###################################
BW   <- 200          # mmol/m3   bottom water conc
Ds   <- 1            # cm2/d     diffusion coeff

O2conc    <- function(Ds,k,BW,x)   BW*exp(-sqrt(k/Ds)*x)
O2flux    <- function(Ds,k,BW)     Ds*sqrt(k/Ds)*BW
O2depth   <- function(Ds,k,BW)     log(BW/0.1)/sqrt(k/Ds)

#======================================
# vertical profiles, different k-values
#======================================

x   <- seq(0,2,length=200)
 
plot (O2conc(Ds,k=10 ,BW,x),x,lty=1,ylim=rev(range(x)),
      xlab="µmol/l",ylab="cm",main="O2",type="l",lwd=2)

lines(O2conc(Ds,k=1  ,BW,x),x,lty=2,lwd=2)
lines(O2conc(Ds,k=100,BW,x),x,lty=3,lwd=2)

legend("right",c("k=10","k=1","k=100"),lty=1:3,lwd=2)
writelabel("A")

# plot flux and penetration depth to screen:
O2flux(Ds,c(1,10,100),BW)
O2depth(Ds,c(1,10,100),BW)

#======================================
# Flux versus oxygen penetration depth
#======================================

# sequence of k-values
kseq    <- exp(seq(log(1),log(1000),length=200))
flux    <- O2flux(Ds,kseq,BW)/100
OPD     <- O2depth(Ds,kseq,BW)*10 
 
plot(flux,OPD,type="l",ylab="Oxygen Penetration depth,mm",xlab="oxygen flux,mmol/m2/d",main="sediment oxygen",lwd=2)
writelabel("B")
#======================================
# calibration
#======================================

# input data
depth <- c(0,0.1,0.2,0.3,0.4,0.5,1.0)
O2    <- c(300,65.52,14.31,3.13,0.68,0.59,0.27)

profiledepth <- c(0,     0.1, 0.2,   0.3,  0.4,   0.5, 1)
O2measured   <- c(300, 65.52, 14.31, 3.13, 0.68, 0.59, 0.27)

BW   <- 300
Ds   <- 1         

plot( O2,depth,lty=1,ylim=c(1,0),pch=15,cex=2,
      xlab="µmol/l",ylab="cm",main="O2")
writelabel("C")
fit <- nls(O2 ~ O2conc(Ds,k,BW,depth),start=c(k=100))   
k   <- coef(fit)

x   <-seq(from=0,to=1,length.out=100)   
lines(O2conc(Ds,k,BW,x),x)

k
O2flux(Ds,k,BW)/100 
