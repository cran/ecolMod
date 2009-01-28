###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 6.                                ##
## Model solution - numerical methods        ##
## Chapter 6.6.5                             ##
## Fate of zooplankton in an estuary         ##
###############################################

#----------------------------------------------------------#
# estuarine advective-diffusive transport and growth/decay #
#----------------------------------------------------------#

require(deSolve)

Zootran <-function(t,Zoo,pars)
{

  with (as.list(pars),{

    Flow   <- meanFlow+ampFlow*sin(2*pi*t/365+phaseFlow)
    seaZoo <- approx(fZooTime, fZooConc, xout=t)$y
    Input  <- +Flow * c(riverZoo, Zoo) +
              -Estar* diff(c(riverZoo, Zoo, seaZoo))
    dZoo   <- -diff(Input)/Volume + g *Zoo
    list(dZoo)
                       })
} 

#-----------------------#
# the forcing function: #
#-----------------------#

# Time and measured value of zooplankton concentration at sea boundary
fZooTime = c(0, 30,60,90,120,150,180,210,240,270,300,340,367)
fZooConc = c(20,25,30,70,150,110, 30, 60, 50, 30, 10, 20, 20)


#-----------------------#
# the model parameters: #
#-----------------------#

# the model parameters:  
pars   <-   c(riverZoo  = 0.0,          # river zooplankton conc
              g         =-0.05,         # /day  growth rate
              meanFlow  = 100*3600*24,  # m3/d, mean river flow
              ampFlow   = 50*3600*24,   # m3/d, amplitude
              phaseFlow = 1.4)          # -     phase of river flow

#--------------------------#
# Initialising morphology: #
#--------------------------#

# parameters defining the morphology
# cross sectional surface area is a sigmoid function of estuarine distance
nbox    <- 100                          
Length  <- 100000                           # m 

dx      <- Length/nbox                      # m

IntDist <- seq(0,by=dx,length.out=nbox+1)   # m
Dist    <- seq(dx/2,by=dx,length.out=nbox)  # m

IntArea <- 4000 + 76000 * IntDist^5 /(IntDist^5+50000^5)   # m2
Area    <- 4000 + 76000 * Dist^5    /(Dist^5+50000^5)      # m2

Volume  <- Area*dx                          # m3


#--------------------------#
# Transport coefficients:  #
#--------------------------#

# parameters defining the dispersion coefficients
# a linear function of estuarine distance

Eriver   <- 0                               # m2/d 
Esea     <- 350*3600*24                     # m2/d 
E        <- Eriver + IntDist/Length * Esea  # m2/d 

Estar  <- E * IntArea/dx                   # m3/d

#----------------------#
# RUNNING the model:   #
#----------------------#
ZOOP  <- rep(0,times=nbox)
times <- 1:365
out   <- ode.band(times=times,y=ZOOP,func=Zootran,parms=pars,nspec=1)

#------------------------#
# PLOTTING model output: #
#------------------------#

# Plot zooplankton; first get the data
par(oma=c(0,0,3,0))   # set margin size (oma)
par(oma=c(0,0,3,0))       # set margin size

filled.contour(x=times,y=Dist/1000,z=out[,-1],
               color= terrain.colors,xlab="time, days",
               ylab= "Distance, km",main="Zooplankton, mg/m3")
mtext(outer=TRUE,side=3,"Marine Zooplankton in the Scheldt",cex=1.5) 
