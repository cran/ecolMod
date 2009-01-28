###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 5.                                ##
## Model solution - analytical methods       ##
## Chapter 5.5. Project                      ##
## Organic matter sinking through a water    ##
## column                                    ##
###############################################

require(shape)

par(mfrow=c(2,2))
par(mar=c(4.1,4.1,5.1,2.1))
Flux  <- 100
v     <- 50
k     <- 0.2
depth <- 0:400


curve(exp(-k/x*400),from=1,to=1000,
      ylab="part deposited at 400m",
      xlab="sinking rate, m/d",lwd=2)
writelabel("A")

curve(exp(-x/v*400),from=0,to=1,
      ylab="part deposited at 400m",
      xlab="decay rate, /d",lwd=2)
writelabel("B")

Conc  <- Flux/v * exp(-k/v*depth)
plot(Conc,depth,ylim=c(400,0),xaxt="n",
     ylab="depth, m",xlab= "",type="l",lwd=2)

axis(side=3)
mtext(side=3,"conc, mmol/m3",line=2,cex=0.8)
writelabel("C")

mtext(side=3,outer=TRUE,"Organic matter sinking through water column",line=-2,cex=1.5)