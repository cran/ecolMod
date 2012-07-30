###############################################
## An NPZD model (Nutrient, Phytoplankton,   ##
## Zooplankton, Detritus)                    ##
##                                           ##
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 2. Model formulation              ##
## 2.9.1                                     ##
###############################################

# load package with the integration routine:

library(deSolve)

#----------------------#
# the model equations: #
#----------------------#

NPZD<-function(t,state,parameters)
 {
 with(as.list(c(state,parameters)),{  # unpack the state variables, parameters
    # Light, a sine function; 50% of light is PAR
    PAR            <- 0.5*(540+440*sin(2*pi*t/365-1.4))
    din            <- max(0,DIN)      # to avoid errors when DIN becomes slightly negative..
    Nuptake        <- maxUptake * PAR/(PAR+ksPAR) * din/(din+ksDIN)*PHYTO

    Grazing        <- maxGrazing* PHYTO/(PHYTO+ksGrazing)*ZOO
    Faeces         <- pFaeces * Grazing
    Excretion      <- excretionRate * ZOO
    Mortality      <- mortalityRate * ZOO * ZOO
    Mineralisation <- mineralisationRate * DETRITUS
    Chlorophyll    <- chlNratio * PHYTO
    TotalN         <- PHYTO + ZOO + DETRITUS + DIN   
  
    dPHYTO    <- Nuptake - Grazing
    dZOO      <- Grazing - Faeces - Excretion - Mortality
    dDETRITUS <- Mortality - Mineralisation + Faeces
    dDIN      <- Mineralisation + Excretion - Nuptake
    
    # the output, packed as a list
    list(c(dPHYTO,dZOO,dDETRITUS,dDIN),                          # the rate of change
        c(Chlorophyll = Chlorophyll, PAR=PAR, TotalN= TotalN))   # the ordinary output variables
    })
  }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#

parameters<-c(maxUptake          =1.0,       # /day
              ksPAR              =140,       # µEinst/m2/s
              ksDIN              =0.5,       # mmolN/m3
              maxGrazing         =1.0,       # /day
              ksGrazing          =1.0,       # mmolN/m3
              pFaeces            =0.3,       # -
              excretionRate      =0.1,       # /day
              mortalityRate      =0.4,       # /(mmolN/m3)/day 
              mineralisationRate =0.1,       # /day
              chlNratio          =1)         # mgChl/mmolN

#-------------------------#
# the initial conditions: #
#-------------------------#
 
state     <-c(PHYTO   =1,                    # state variable initial conditions, units mmolN/m3
              ZOO     =0.1,
              DETRITUS=5.0,
              DIN     =5.0)

#----------------------#
# RUNNING the model:   #
#----------------------#
#  2 steps
#  step 1  : run the model for 365 days, no intermediate output required

times     <-c(0,365)                        
out       <-as.data.frame(ode(state,times,NPZD,parameters)) 

#  step 2  : the model is reinitialised with the final conditions of previous run

num       <- length(out$PHYTO)   # last element 
state     <- c(PHYTO   =out$PHYTO[num],         
               ZOO     =out$ZOO[num],           
               DETRITUS=out$DETRITUS[num],
               DIN     =out$DIN[num])

# new run; this time the output is kept: time: from 1 to 365 days, outputsteps of days
times     <-seq(0,730,by=1)              
out       <-as.data.frame(ode(state,times,NPZD,parameters))

#------------------------#
# PLOTTING model output: #
#------------------------#
windows()
par(mfrow=c(2,2), oma=c(0,0,3,0))   # set number of plots (mfrow) and margin size (oma)

plot (times,out$PAR        ,type="l",lwd=2,main="PAR"        ,xlab="time, hours",ylab="µEinst/m2/s")
#writelabel("A") needs package (shape)
plot (times,out$Chlorophyll,type="l",lwd=2,main="Chlorophyll",xlab="time, hours",ylab="µg/l")
#writelabel("B")
plot (times,out$ZOO        ,type="l",lwd=2,main="Zooplankton",xlab="time, hours",ylab="mmolN/m3")
#writelabel("C")
plot (times,out$DIN        ,type="l",lwd=2,main="DIN"        ,xlab="time, hours",ylab="mmolN/m3")
#writelabel("D")

mtext(outer=TRUE,side=3,"NPZD model",cex=1.5)
