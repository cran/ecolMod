###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 3.                                ##
## Spatial components and transport          ##
## Chapter 3.6.2.                            ##
## A 1-D microscopic and macroscopic model   ##
## of diffusion.                             ##
###############################################

windows()
par(mfcol=c(2,2))
require (shape)

#======================
# microscopic model
#======================

nind    <- 100 
nsteps  <- 100
pos     <- matrix (nrow=nsteps+1,ncol=nind,data=0)

# random steps, between -0.5 and 0.5
for (i in 1:nsteps) pos[i+1,] <- pos[i,]+runif(nind)-0.5


matplot(pos,type="l",col="black",lty=1,xlab="step",ylab="Position")
mtext("microscopic model",side=3,line=2,cex=1.2)

# density plot
posmicro <- matrix(ncol=nsteps+1,nrow=100)
for (i in 1:(nsteps+1)) 
  posmicro[,i] <- density(pos[i,],from=-10,to=10,n=100)$y


yy     <- 1:nsteps
xx     <- seq(-10,10,length.out=100)

matplot(posmicro[,seq(1,100,by=4)],type="l",lty=1,col="black",
        axes=FALSE,frame.plot=TRUE,ylab="Density")

#=============================
# macroscopic model
#=============================
# estimate the parameters of the continuous time model:
vari <- apply(X=pos,MARGIN=1,FUN=var)
l1   <- lm(vari~ c(1:101) +0)

Ds  <-  coef(l1)/2   #0.075/2   # diffusion coefficient
ini <- 1       # initial condition


# analytical solution of the 1-D diffusion equation in cartesian coordinates  
xx <-seq(-10,10,length=100)
tt <-seq(1,100,by=1)
posmacro <- outer(xx,tt,
    FUN =function (xx,tt) ini/(2*sqrt(pi*Ds*tt))*exp(-xx^2/(4*Ds*tt)))

persp(xx,tt,z=posmacro,theta=150,box=TRUE,axes=TRUE,border=NA,
     xlab="position",ylab="time",zlab="density",col=drapecol(posmacro))  
mtext("macroscopic model",side=3,line=2,cex=1.2)

matplot(posmacro[ ,seq(1,100,by=4)],type="l",lty=1,col="black",
        axes=FALSE,frame.plot=TRUE,ylab="Density")

