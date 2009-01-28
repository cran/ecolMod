###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.8.5                             ##
## Fate of marine zooplankton in an estuary  ##
## equilibrium condition                     ##
###############################################

# load package with the steady-state routines:
require(rootSolve)

Estuary <-function(t,C,pars)
{

  with (as.list(pars),{

    Input <-  Q     * c(Criver, C) +
             -Estar * diff(c(Criver, C, Csea))
    dC    <- -diff(Input)/Volume + rate*C
    list(dC)
  })
}


nbox    <- 100
Length  <- 100000                       # m total length of estuary
dx      <- Length/nbox

# distance from river to interfaces and to mid of boxes  (m)
IntDist   <- seq(0,by=dx,length.out=nbox+1)
Dist      <- seq(dx/2,by=dx,length.out=nbox)

# Sigmoidal increase in cross sectional area
# Cross section at interfaces and midle of boxes         (m2)
Area_int <- 4000 + 76000 * IntDist^5 /(IntDist^5+50000^5)
Area     <- 4000 + 76000 * Dist^5    /(Dist^5+50000^5)

# Volume of boxes                                        (m3)
Volume   <- Area*dx

Eriver   <- 0                      # m2/d dispersion coefficient at upstream boundary
Esea     <- 350*3600*24            # m2/d dispersion coefficient at seaward boundary
E        <- Eriver + IntDist/Length * Esea

# BulkDisp is used in the calculations; m3/d
Estar  <- E * Area_int/dx

# the model parameters:
pars   <-   c(Criver  = 0.0,          # riverine conc
              Csea    = 100.0,        # seaward conc
              rate    = -0.05,        # /day  growth rate
              Q       = 100*3600*24)  # m3/d, mean river flow

pars["rate"] <- 0
pars["Csea"] <- 35

sal <- steady.band(y=rep(35,times=nbox),func=Estuary,
                   parms=pars,nspec=1)$y

Zoo  <- NULL

gSeq <-seq (-0.05,0.01, by=0.01)
for (g in gSeq)
{
pars["rate"] <- g
pars["Csea"] <- 100
st <- steady.band(y=rep(100,nbox),func=Estuary,parms=pars,nspec=1 ,pos=TRUE  ) #
Zoo <-  cbind(Zoo,st$y)
}

matplot(sal,Zoo,type="l",lwd=1,col="black",xlab="Salinity",
        ylab="gDWT/m3",main="Zooplankton",lty=1:20)
legend("topleft",legend= gSeq,title="g, /day",lty=1:20)

