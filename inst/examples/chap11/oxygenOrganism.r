###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 11.                               ##
## Testing and validating the model          ##
## Chapters 11.6.s                           ##
## Time-varying Oxygen consumption in a      ##
## small cylindrical organism                ##
###############################################

require (deSolve)
require (rootSolve)
require (shape)

# parameters

BW     <- 2            # mmol/m3,       oxygen concentration in surrounding water
Da     <- 0.5          # cm2/d          effective diffusion coefficient in organism
R      <- 0.0025        # cm             radius of organism
Q      <- 250000       # nM/cm3/d       oxygen consumption rate per volume per day
L      <- 0.05          # cm             length of organism

# the analytical cylindrical model
cylinder <- function(Da,Q,BW,R,r)  BW+Q/(4*Da)*(r^2-R^2)

par(mfrow=c(2,2))

N  <- 10                           # we consider 10 layers in the nematode body
dx <- R/N                          # thickness of each layer
x  <- seq(dx/2,by=dx,length.out=N) # distance of center to mid-layer
xi <- seq(0,by=dx,length.out=N+1)  # distance to layer interface
dxi<- c(rep(dx,N),dx/2)            # dispersion distances
A   <- 2*pi*x *L                   # surface at mid-layer depth
Ai  <- 2*pi*xi*L                   # surface at layer interface

# the numerical model

oxygen <- function (time, O2, pars)
#==============================================================================
# the rate of change of oxygen in the body of a cylindrical organism
# outer boundary = concentration boundary
# lower boundary = zero-gradient boundary
#==============================================================================

  {
    BWO2 <- BW*(1-0.8*sin(2*pi*time*24))
    # outer concentration imposed (BW), lower: zero gradient
    Flux    <- -Da * diff(c(O2[1],O2,BWO2))/dxi   #diffusion

    # Rate of change = Flux gradient - oxygen consumption
    dO2     <- -diff(Ai * Flux)/A/dx -Q

    return (list(dO2=dO2,c(Flux=Flux,BWO2=BWO2)))
  }


CONC  <- steady.1D (runif(N),func=oxygen,nspec=1,atol=1e-10)
O2    <- CONC$y
plot(x,O2,xlab="distance, cm",ylab="oxygen, µmol/l")
lines(x, BW+Q/(4*Da)*(x^2-R^2))
legend ("topleft",lty=c(1,NA),pch=c(NA,1),c("analytical solution","numerical approximation"))
writelabel("A")

N  <- 100                          # we consider 100 layers in the nematode body
dx <- R/N                          # thickness of each layer
x  <- seq(dx/2,by=dx,length.out=N) # distance of center to mid-layer
xi <- seq(0,by=dx,length.out=N+1)  # distance to layer interface
dxi<- c(rep(dx,N),dx/2)            # dispersion distances
A   <- 2*pi*x *L                   # surface at mid-layer depth
Ai  <- 2*pi*xi*L                   # surface at layer interface

CONC  <- steady.1D (runif(N),func=oxygen,nspec=1,atol=1e-10)
O2    <- CONC$y
plot(x,O2,xlab="distance, cm",ylab="oxygen, µmol/l")
lines(x, BW+Q/(4*Da)*(x^2-R^2))
legend ("topleft",lty=c(1,NA),pch=c(NA,1),c("analytical solution","numerical approximation"))
writelabel("B")


times<- seq(0,1/24,length.out=120)
out   <-as.data.frame(ode.1D(O2,times,func=oxygen,parms=0,rtol=1e-10,atol=1e-10,nspec=1))
oxy   <- out[,2:101]
plot(times*24,out$BWO2,xlab="time, hour",ylab="µmol/l", main="BW concentration",type="l", lwd=2)
writelabel("C")
image(times*24,y=x,z=as.matrix(oxy),xlab="time, hour",ylab="distance ",main="Dynamic simulation",col=grey(seq(0.1,1,len=100)))
contour(times*24,y=x,z=as.matrix(oxy),add=TRUE)
writelabel("D")
mtext(outer=TRUE,side=3,"Oxygen in cylindrical body",cex=1.5,line=-1.5)

