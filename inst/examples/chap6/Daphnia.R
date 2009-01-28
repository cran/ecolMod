###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 6.                                ##
## Model solution - numerical methods        ##
## Chapter 6.6.2                             ##
## Growth of a Daphnia individual            ##
###############################################

# load package with the integration routine:

library(deSolve)

#----------------------#
# the model equations: #
#----------------------#

model<-function(t,state,parameters)
 {
 with(as.list(c(state)),{  # unpack the state variables

  # ingestion, size-dependent and food limited
  WeightFactor <- (IngestWeight-INDWEIGHT)/(IngestWeight-neonateWeight) 
  MaxIngestion <- maxIngest*WeightFactor      # /day     
  Ingestion    <- MaxIngestion*INDWEIGHT*FOOD / (FOOD + ksFood)

  Respiration  <- respirationRate * INDWEIGHT         # µgC/day 
  Growth       <- Ingestion*assimilEff - Respiration

  # Fraction of assimilate allocated to reproduction

  if (Growth <= 0. | INDWEIGHT<reproductiveWeight) Reproduction <- 0. 
  else {               # Fraction of growth allocated to reproduction.
    WeightRatio  <- reproductiveWeight/INDWEIGHT
    Reproduction <- maxReproduction * (1. - WeightRatio^2)
       }

  # rate of change 
  dINDWEIGHT <- (1. -Reproduction) * Growth
  dEGGWEIGHT <-      Reproduction  * Growth
  dFOOD      <- -Ingestion * numberIndividuals      

  # the output, packed as a list
    list(c(dINDWEIGHT, dEGGWEIGHT, dFOOD),   # the rate of change
         c(Ingestion    = Ingestion,            # the ordinary output variables
           Respiration  = Respiration,
           Reproduction = Reproduction))
    })

  }  # end of model

#----------------------#
# Moulting weight loss #
#----------------------#

Moulting   <- function ()

  {
   with(as.list(c(state)),{  # unpack the state variables
 
   # Relationship moulting loss and length
    refLoss   <-  0.24   #µgC
    cLoss     <-  3.1    #-

    # Weight lost during molts depends allometrically on the organism length
    INDLength    <- (INDWEIGHT /3.0)^(1/2.6)

    WeightLoss <- refLoss * INDLength^cLoss
    return(INDWEIGHT - WeightLoss)   # New weight
    })
  }

#-----------------------#
# the model parameters: #
#-----------------------#

neonateWeight      <-  1.1    #µgC
reproductiveWeight <-  7.5    #µgC
maximumWeight      <- 60.0    #µgC

ksFood             <- 85.0    #µgC/l
IngestWeight       <-132.0    #µgC
maxIngest          <-  1.05   #/day
assimilEff         <-  0.8    #-

maxReproduction    <-  0.8    #-
respirationRate    <-  0.25   #/day

# Dilution parameters !
transferTime       <-    2    # Days
foodInMedium       <-  509    # µgC/l

instarDuration     <-  3.0    # days
numberIndividuals  <-   32    #   -

#-------------------------#
# the initial conditions: #
#-------------------------#
 
state     <-c(
  INDWEIGHT = neonateWeight      , # µgC
  EGGWEIGHT = 0                  , # µgC    ! Total egg mass in a stage
  FOOD      = foodInMedium         # µgC
             )

#----------------------#
# RUNNING the model:   #
#----------------------#

TimeFrom     <- 0
TimeEnd      <- 40                          # duration of simulation, days
TimeMoult    <- TimeFrom + instarDuration   # next time (days) at which moulting 
TimeTransfer <- TimeFrom + transferTime     # next time (days) at which individuals are transferred

Time         <- TimeFrom
Outdt        <- 0.1                         # output time step
out          <- NULL                        # output array

while (Time < TimeEnd)
{
  TimeOut <- min(TimeMoult,TimeTransfer,TimeEnd)  # integrator runs till TimeOut
  times   <- seq(Time,TimeOut,by=Outdt)           # sequence of output times
  if (length(times)>1) {
  out1    <-as.data.frame(ode(state,times,model,parms=0))  # integrate
  out     <- rbind(out,out1)                                 # add output to output array
  lout    <- nrow(out1)                                      # last element of output
  state   <-c(
             INDWEIGHT = out1[lout,"INDWEIGHT"],
             EGGWEIGHT = out1[lout,"EGGWEIGHT"],
             FOOD      = out1[lout,"FOOD"])
  }
  if (Time >= TimeMoult)     # Moulting...
    {
     state[1]     <- Moulting()  # New weight individuals
     state[2]     <- 0.          # Put eggs = 0
     TimeMoult    <- Time +instarDuration          # next time at which moulting occurs
    }
  if (Time >= TimeTransfer)  # New medium...
    {
     state[3]     <- foodInMedium
     TimeTransfer <- Time + transferTime           # next time at which individuals are transferred
    }
  
  # Reset time, state variables
  Time   <- TimeOut 

 }


#------------------------#
# PLOTTING model output: #
#------------------------#

windows()
par(mfrow=c(2,2), oma=c(0,0,3,0))   # set number of plots (mfrow) and margin size (oma)

plot (out$time,out$FOOD        ,type="l",main="Food"              ,xlab="time, days",ylab="gC/m3")
plot (out$time,out$INDWEIGHT   ,type="l",main="individual weight" ,xlab="time, days",ylab="µgC")
plot (out$time,out$EGGWEIGHT   ,type="l",main="egg weight"        ,xlab="time, days",ylab="µgC")
plot (out$time,out$Ingestion   ,type="l",main="Ingestion"             ,xlab="time, days",ylab="µgC/day")

mtext(outer=TRUE,side=3,"DAPHNIA model",cex=1.5)

windows()
par (mfrow=c(2,2))
curve(maxIngest*(IngestWeight-x)/(IngestWeight-neonateWeight),0,60,
      main="Max. ingestion rate",ylab="/d",xlab="ind. weight, µC",lwd=2) 
curve(pmax(0., maxReproduction * (1. - (reproductiveWeight/x)^2)),0,60,
      main="fraction assimilate to reproduction ",ylab="-",
      xlab="ind. weight, µC",lwd=2) 
curve(((x /3.0))^(1/2.6),0,60,
      main="Individual length",ylab="µm",xlab="ind. weight, µC",lwd=2) 
curve(0.24*((x /3.0)^(1/2.6))^3.1,0,60,
      main="Weight loss during moulting",ylab="µg",xlab="ind. weight, µC",lwd=2) 


