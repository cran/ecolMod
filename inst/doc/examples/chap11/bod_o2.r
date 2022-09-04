###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 11.                               ##
## Testing and validating the model          ##
## Chapters 11.3-11.5                        ##
## The BOD-Oxygen model                      ##
###############################################

# load package with the integration routine:
require(deSolve)
require(shape)
par(mfrow=c(2,2))

k       = 0.1             # /day     - reaeration
O2sat   = 300             # mmol/m3  - saturated oxygen concentration
r       = 0.05            # /day     - BOD decay rate
O2_0    = 250             # mmol/m3  - Initial oxygen concentration
BOD_0   = 500             # mmol/m3  - Initial BOD concentration
ks      = 0               # mmol/m3  - half-saturation concentration

# numerical model
numBOD <- function (time,state,pars)
{
 with (as.list(state),
  {
    dO2  <- -r*BOD*O2/(O2+ks)+k*(O2sat-O2)
    dBOD <- -r*BOD*O2/(O2+ks)
    return(list(c(dO2,dBOD)))
  }
      )
}

# A comparison numerical / analytical model
# numerical solution plotted as points
times <- 0:100
state <- c(O2=O2_0,BOD=BOD_0)
out   <- as.data.frame(ode(state,times,numBOD,0))
plot(out$time,out$O2,xlab="time",ylab="mmol O2/m3",lwd=2,main="Correctness of solution") 

# analytical solution - added as a curve
curve(analytical(x,k),lty=1,lwd=1,add=TRUE)
legend("bottomright",c("analytic","numerical"),lwd=c(1,2),lty=c(1,NA),pch=c(NA,1))
writelabel("A")

# B: internal logic
# wrong use of model : too low reaeration -> negative concentration

k     <- 0.01
times <- 0:200
state <- c(O2=O2_0,BOD=BOD_0)
out   <- as.data.frame(ode(state,times,numBOD,0))
plot(out$time,out$O2,xlab="time",ylab="mmol O2/m3",main="Internal logic",type="l",lty=2) 
abline(h=0,lty=3)

ks    <- 1
state <- c(O2=O2_0,BOD=BOD_0)
out2  <- as.data.frame(ode(state,times,numBOD,0))
lines(out2$time,out2$O2,lwd=2)
legend("bottomright",c("no O2 limitation","O2 limitation"),lwd=c(1,2),lty=c(2,1))
writelabel("B")

# C: global sensitivity
k       <- 0.1           
rseq   <- seq(0.0,0.2,by=0.002)
rseq   <- rseq[rseq!=k]    # cannot calculate analytical solution for this...
minO2  <- rseq
for (i in 1:length(rseq)) 
  minO2[i]  <- min(analytical(times,r=rseq[i]))
plot(rseq,minO2,type="l",lwd=2 ,xlab="r, /day",ylab="minimum O2, mmol/m3",main="global sensitivity")
writelabel("C")

mtext(side=3,outer=TRUE,line=0,"BOD-O2 model",cex=1.25,font=2)

# D: local sensitivity

times   <- 0:100 
ss      <- 1.1

kp      <- k * ss         # /day     - reaeration
O2satp  <- O2sat*ss        # mmol/m3  - saturated oxygen concentration
rp      <- r*ss           # /day     - BOD decay rate

ref  <- analytical(times)
outk <- analytical(times,k=kp)
outs <- analytical(times,O2sat=O2satp)
outr <- analytical(times,r=rp)
outm <- mean(ref)
ss   <- cbind(k=(outk-ref)/outm/0.1,sat=(outs-ref)/outm/0.1,r=(outr-ref)/outm/0.1)

plot(times,ref,ylim=range(c(ref,outs)),type="l",lwd=2,xlab="time",ylab="mmol O2/m3",main="local sensitivity")
lines(times,outs,lwd=2,lty=2)
arrseq <- seq(10,100,10)#c(10,30,50,70,90)
Arrows(times[arrseq],ref[arrseq],times[arrseq],outs[arrseq],arr.len=0.25,arr.adj=1)
legend("topleft",c(expression(O[2]^"*"== 300),expression(O[2]^"*"== 330)),lwd=2,lty=c(1,2))

writelabel("D")
