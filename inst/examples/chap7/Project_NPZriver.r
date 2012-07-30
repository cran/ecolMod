###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.9.5. Project                    ##
## Succession of nutrients, phytoplankton,   ##
## and zooplankton in a river                ##
###############################################
require(rootSolve)
#==============================================================================
# NPZmodel equations
#==============================================================================

NPZmodel <- function (time=0,state,pars=NULL)
{
 N<- state[1       :Nb    ]
 P<- state[(Nb+1)  :(2*Nb)]
 Z<- state[(2*Nb+1):(3*Nb)]

# transport          #
# advective fluxes at upper interface of each layer
# freshwater concentration imposed 

 NFlux <- flow * c(RiverN ,N) 
 PFlux <- flow * c(RiverP ,P) 
 ZFlux <- flow * c(RiverZ ,Z) 
 
# Biology            #

 Pprod <-  mumax * N/(N+kN)*P      # primary production
 Graz  <-  gmax * P/(P+kP)*Z       # zooplankton grazing
 Zmort <-  mrt*Z                   # zooplankton mortality

# Rate of change= Flux gradient and biology
 dN        <- -diff(NFlux)/delx - Pprod + Graz*(1-eff) + Zmort
 dP        <- -diff(PFlux)/delx + Pprod - Graz 
 dZ        <- -diff(ZFlux)/delx         + Graz*eff     - Zmort 

 return(list(c(dN=dN,dP=dP,dZ=dZ),
             c(Pprod=Pprod,Graz=Graz,Zmort=Zmort,
            Nefflux=NFlux[Nb],Pefflux=PFlux[Nb],Zefflux=ZFlux[Nb])))
}


###############################################################################
# model application
###############################################################################
# model parameters
riverlen <- 100               # total length river,         km
Nb       <- 100               # number boxes
delx     <- riverlen / Nb     # box length

RiverN   <- 100               # N concentration at river,   mmolN/m3
RiverP   <- 10                # P concentration at river,   mmolN/m3
RiverZ   <- 1                 # Z concentration at river,   mmolN/m3

flow     <- 1                  # river flow,                 km/day

mumax    <- 0.5               # maximal light-limited primary production, /day
kN       <- 1                 # half-saturated N for pprod, mmolN/m3
gmax     <- 0.5               # max grazing rate,           /day
kP       <- 1                 # half-saturated P for graz , mmolN/m3
mrt      <- 0.05               # mortality rate,             /day
eff      <- 0.7               # growth efficiency,          -

# initial guess of NPZ; just random numbers
Conc     <- runif(3*Nb) #c(N,P,Z)

# steady-state solution, v= 1 km/day
sol  <- steady.1D (y=Conc, func=NPZmodel,nspec=3,maxiter=100,atol=1e-10,positive=TRUE) 
Conc <- sol$y
N    <- Conc[1       :Nb    ]
P    <- Conc[(Nb+1)  :(2*Nb)]
Z    <- Conc[(2*Nb+1):(3*Nb)]

Res  <- NPZmodel(state=Conc)

# run with v=5 km d-1
flow     <- 5                  # river flow,                 km/day
Conc     <- runif(3*Nb) #c(N,P,Z)
sol2  <- steady.1D (y=Conc, func=NPZmodel,nspec=3,maxiter=100,atol=1e-10,positive=TRUE) 
Conc2 <- sol2$y
N2    <- Conc2[1       :Nb    ]
P2    <- Conc2[(Nb+1)  :(2*Nb)]
Z2    <- Conc2[(2*Nb+1):(3*Nb)]

# run with v=10 km d-1
flow     <- 10                  # river flow,                 km/day
Conc     <- runif(3*Nb) #c(N,P,Z)
sol3  <- steady.1D (y=Conc, func=NPZmodel,nspec=3,maxiter=100,atol=1e-10,positive=TRUE) 
Conc3 <- sol3$y
N3    <- Conc3[1       :Nb    ]
P3    <- Conc3[(Nb+1)  :(2*Nb)]
Z3    <- Conc3[(2*Nb+1):(3*Nb)]

windows()
par(mfrow=c(2,2))
Dist <- seq(delx/2,riverlen,by=delx)
plot(Dist,N,ylab="mmolN/m3",main="DIN",type="l",lwd=2,ylim=c(0,110))
lines(Dist,N2)
lines(Dist,N3,lty=2)
plot(Dist,P,ylab="mmolN/m3" ,main="Phytoplankton",type="l",lwd=2,ylim=c(0,110))
lines(Dist,P2)
lines(Dist,P3,lty=2)
plot(Dist,Z,ylab="mmolN/m3" ,main="Zooplankton",type="l",lwd=2,ylim=c(0,110))
lines(Dist,Z2)
lines(Dist,Z3,lty=2)

plot(0,type="n",xlab="",ylab="",axes=FALSE)
legend("center",lty=c(1,1,2),lwd=c(2,1,1),
c(expression(v==1~km~d^{-1}),expression(v==5~km~d^{-1}),expression(v==10~km~d^{-1}))) 

mtext(side=3,outer=TRUE,"NPZ model",line=-1.5,cex=1.5)