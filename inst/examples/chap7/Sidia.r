###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.8.4                             ##
## Steady-state of the Si diagenetic model   ##
###############################################

#======================
# Silicate diagenesis
#======================

require(rootSolve)

SiDiamodel <- function (time=0,Conc,pars=NULL)
{
 BSi<- Conc[1:N]
 DSi<- Conc[(N+1):(2*N)]

  # diffusive fluxes at upper interface of each layer
  # upper concentration imposed (bwDSi), lower: zero gradient
 DSiFlux <- -SedDisp *   IntPor *diff(c(bwDSi ,DSi,DSi[N]))/thick    
 BSiFlux <- -Db      *(1-IntPor)*diff(c(BSi[1],BSi,BSi[N]))/thick 

 BSiFlux[1] <- BSidepo                         # upper boundary flux is imposed

 # BSi dissolution    #

 Dissolution <- rDissSi * BSi*(1.- DSi/EquilSi )^pow 
 Dissolution <- pmax(0,Dissolution)

 # Rate of change= Flux gradient, corrected for porosity and dissolution
 dDSi        <- -diff(DSiFlux)/thick/Porosity      +           # transport
                 Dissolution * (1-Porosity)/Porosity           # biogeochemistry

 dBSi        <- -diff(BSiFlux)/thick/(1-Porosity)  - Dissolution				

 return(list(c(dBSi=dBSi,dDSi=dDSi),
             Dissolution=Dissolution,
             DSiSurfFlux =DSiFlux[1],DSIDeepFlux =DSiFlux[N+1],
             BSiDeepFlux =BSiFlux[N+1]))
}

# sediment parameters
thick    <- 0.05                       # thickness of sediment layers (cm)
Intdepth <- seq(0,10,by=thick)         # depth at upper interface of each layer 
Nint     <- length(Intdepth)           # number of interfaces
Depth    <- 0.5*(Intdepth[-Nint] +Intdepth[-1]) # depth at middle of each layer
N        <- length(Depth)                       # number of layers

por0    <- 0.9                         # surface porosity (-)
pordeep <- 0.7                         # deep porosity    (-)
porcoef <- 2                           # porosity decay coefficient  (/cm)
Porosity <- pordeep + (por0-pordeep)*exp(-Depth*porcoef)     # porosity profile, middle of layers
IntPor   <- pordeep + (por0-pordeep)*exp(-Intdepth*porcoef)  # porosity profile, upper interface

dB0      <- 1/365           # cm2/day       - bioturbation coefficient
dBcoeff  <- 2               
mixdepth <- 5                # cm
Db       <- pmin(dB0,dB0*exp(-(Intdepth-mixdepth)*dBcoeff))

# biogeochemical parameters
SedDisp  <- 0.4              # diffusion coefficient, cm2/d  
rDissSi  <- 0.005            # dissolution rate, /day
EquilSi  <- 800             # equilibrium concentration
pow      <- 1
BSidepo  <- 0.2*100          # nmol/cm2/day
bwDSi    <- 150              # µmol/l  

# initial guess of state variables-just random numbers between 0,1
Conc     <- runif(2*N) 

# three runs with different deposition rates
BSidepo  <- 0.2*100          # nmol/cm2/day
sol  <- steady.1D (y=Conc, func=SiDiamodel, parms=NULL, nspec=2) 
CONC <- sol$y

BSidepo  <- 2*100          # nmol/cm2/day
sol2 <- steady.1D (y=Conc, func=SiDiamodel,parms=NULL, nspec=2) 
CONC2 <- sol2$y

BSidepo  <- 3*100          # nmol/cm2/day
sol3 <- steady.1D (y=Conc, func=SiDiamodel,parms=NULL, nspec=2) 
CONC3 <- sol3$y

DSi  <- cbind(CONC[(N+1):(2*N)],CONC2[(N+1):(2*N)],CONC3[(N+1):(2*N)])
BSi  <- cbind(CONC[1:N],CONC2[1:N],CONC3[1:N])

windows()
par(mfrow=c(2,2))

matplot(DSi,Depth,ylim=c(10,0),xlab="mmolSi/m3 Liquid",main="DSi",type="l",lwd=c(1,2,1),col="black")
matplot(BSi,Depth,ylim=c(10,0),xlab="mmolSi/m3 Solid" ,main="BSi",type="l",lwd=c(1,2,1),col="black")
legend("right",c("0.2","2","3"),title="mmol/m2/d",lwd=c(1,2,1),lty=1:3)
plot(Porosity,Depth,ylim=c(10,0),xlab="-" ,main="Porosity",type="l",lwd=2)
plot(Db,Intdepth,ylim=c(10,0),xlab="cm2/d" ,main="Bioturbation",type="l",lwd=2)
