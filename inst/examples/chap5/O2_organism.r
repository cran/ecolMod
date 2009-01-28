###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 5.                                ##
## Model solution - analytical methods       ##
## Chapter 5.4.3.                            ##
## Steady-state oxygen budget in a small     ##
## organism living in suboxic conditions     ##
###############################################

windows()
par(mfrow=c(2,2))
require(shape)

##########################################
# 1-dimensional oxygen model in organism
##########################################

# parameters

BW     <- 2            # mmol/m3,       oxygen concentration in surrounding water
Da     <- 0.5          # cm2/d          effective diffusion coefficient in organism
R      <- 0.005        # cm             radius of organism 
Q      <- 250000       # nM/cm3/d       oxygen consumption rate per volume per day

# distance in organism body
rr        <- seq(-R,R,length=400)

# sandwich model (infinite plate), the cylindrical model and the spherical model
sandwich <- function(Da,Q,BW,R,r)  BW+Q/(2*Da)*(r^2-R^2)
cylinder <- function(Da,Q,BW,R,r)  BW+Q/(4*Da)*(r^2-R^2)
sphere   <- function(Da,Q,BW,R,r)  BW+Q/(6*Da)*(r^2-R^2)

# Oxygen concentration in organisms body, for 3 scenarios; BW conc=2,width=100µm
plot  (rr,sandwich(Da,Q,BW,R,rr),lwd=2,lty=1,
       ylim=c(0,BW), type="l",main="oxygen in organism" ,
       xlab="radius,cm",ylab="µmol/l")
lines (rr,cylinder(Da,Q,BW,R,rr),lwd=2,lty=2)
lines (rr,sphere  (Da,Q,BW,R,rr),lwd=2,lty=3)
legend("top",c("sheet","cylinder","sphere"),lty=1:3,lwd=2)
writelabel("A")

# oxygen concentration for a cylinder, body width varying
plot  (rr,cylinder(Da,Q,BW,R,rr),ylim=c(0,BW), type="l",lwd=2,
       xlab="radius,cm",ylab="O2, µmol/l",main="cylinder" )
lines (rr*0.5,cylinder(Da,Q,BW,R*0.5,rr*0.5),lwd=2,lty=2)
lines (rr*0.25,cylinder(Da,Q,BW,R*0.25,rr*0.25),lwd=2,lty=3)
legend("bottom",legend=c("100µm","50µm","25µm"),title="maximal thickness",lty=1:3,lwd=2)
writelabel("B")

#Functions to estimate the critical thickness
critsand <- function(Da,Q,BW)  sqrt(BW*2*Da/Q)
critcyl  <- function(Da,Q,BW)  sqrt(BW*4*Da/Q)
critsph  <- function(Da,Q,BW)  sqrt(BW*6*Da/Q)

BWseq <- seq(0,50,length=100)
plot(BWseq,10000*sqrt(BWseq*6*Da/Q),type="l",lty=3,lwd=2,
    main="critical thickness", xlab="surrounding oxygen, µmol/l",ylab="µm")
lines(BWseq,10000*sqrt(BWseq*4*Da/Q),lty=2,lwd=2)
lines(BWseq,10000*sqrt(BWseq*2*Da/Q),lty=1,lwd=2)
legend("bottomright",c("sheet","cylinder","sphere"),lty=1:3,lwd=2)

writelabel("C")

emptyplot(c(-1.5,1.5),main="shapes" )
ypos <- c(0.9,1.1)
plotellipse(mid=c(0,ypos[1]),rx=1,ry=0.4,from=-pi,to=0,lwd=1)

plotellipse(mid=c(0,ypos[2]),rx=1,ry=0.4,col="lightgrey",lwd=1)
plotellipse(mid=c(0,ypos[1]),rx=1,ry=0.4,from=0,to=pi,lwd=1,lcol="darkgrey")
segments(-1,ypos[1],-1,ypos[2],lwd=2)
segments( 1,ypos[1], 1,ypos[2],lwd=2)

rseq1 <- seq(0.3,1,length.out=50)
col  <- grey(c(rseq1,rev(rseq1)))
filledcylinder(rx=0.1,ry=0.2,len=2.5,angle=0,col=col,mid=c(0,0),lcol="black",lwd=1)

col2  <- grey(seq(1,0.5,length.out=50))
filledcircle(r1=0.5,col=col2,mid=c(0,-1.),lcol="black",lwd=1)
box()
writelabel ("D")
