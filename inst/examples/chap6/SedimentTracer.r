###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 6.                                ##
## Model solution - numerical methods        ##
## Chapter 6.5.1 Sediment tracer model:      ##
## numerical diffusion                       ##
###############################################

require(deSolve)

# general parameters
times     <- c(1/365,10/365,30/365,1/4,1/2,1)          # y    time at which output is required

# General parameters for all runs
DD        <- 0.1      # cm2/y   diffusion coefficient representing bioturbatin
u         <- 1       # cm/y    advection = sedimentation rate
init      <- 1000    # n.cm-3  initial density in upper layers
outtime   <- 6       # choose the output time point (as found in 'times') you want to see in the graphs
                     # default outtime=6 will show distribution after 1.0 year

# Function to calculate derivatives
Lumin <-function(t,LUMIN,parameters)
 {  with (as.list(c(LUMIN,parameters)),
    {
    FluxDiff     <- -DD*diff(c(LUMIN[1],LUMIN,LUMIN[numboxes]))/delx
    FluxAdv      <- u * (sigma * c(0,LUMIN) + (1-sigma) * c(0,LUMIN[2:numboxes],LUMIN[numboxes]))
    Flux         <- FluxDiff + FluxAdv

    dLUMIN       <- -diff(Flux)/delx
    
    list(dLUMIN )
     })
 }
 
# Function to run the model and produce output

Modelfunc <- function(ltt,diffm,numboxes)
{ 
    # complete set of parameters, depending on choices 
     delx      <- 3/numboxes                                     # thickness of boxes in cm
    if (diffm=="back")   sigma<-1
    if (diffm=="cent")   sigma<-0.5
    if (diffm=="fiad") {
                        if (DD>0) { Pe <- u * delx/DD
                                    sigma <- (1 + (1/tanh(Pe) - 1/Pe))/2
                                  }
                        if (DD==0) sigma <- 1
                       }
     pars<-as.list(c(          
              numboxes = numboxes,
              delx     = delx,
              sigma    = sigma,
              DD       = DD,
              u        = u
              ))

    # initial conditions
     borders         <- seq(from=0,by=delx,length.out=numboxes) 
     Ll              <- rep(0,times=numboxes)       # ind/m2
     for (i in 1:numboxes-1){     
       if (borders[i+1]<=0.5)                       Ll[i]  <- init
       if ((borders[i]<0.5) && (borders[i+1]>0.5))  Ll[i]  <- init * (0.5-borders[i])/delx
                            }

    # call
     state           <- c(LUMIN=Ll)
     out             <- ode.band(y=state,times,func=Lumin,parms=pars,nspec=1)
    # store output
     Dl              <- out[outtime,2:(numboxes+1)]
    # output
     Distance  <- seq(from=delx/2, by=delx, length.out=numboxes) # distance of box centres from x=0
     if (diffm=="back") ttext<-"Backward Differences"
     if (diffm=="cent") ttext<-"Centered Differences"
     if (diffm=="fiad") ttext<-"Fiadeiro Scheme"
     if (ltt==1){ 
        plot(Dl,Distance,xlim=c(0,max(Dl)),ylim=c(numboxes*delx,0),type="l",
           main=ttext,xlab="Tracer density",ylab="Depth (cm)")
        legend(200,2, legend=c("n=300","n=120","n=60","n=30","n=15"),  
              lty=c(1,2,3,4,5),title="no. of boxes",pch=c(-1,-1,-1,-1,-1))
                } 
     if (ltt>1)lines(Dl,Distance,lty=ltt)
}

# Run model with number of different options
    windows()
    par (mfrow=c(1,3),oma=c(1,1,4,1),mar=c(5, 4, 4, 2)+0.1,mgp=c(3,1,0)) 
    numboxlist <- c(300,120,60,30,15)
    
    for (diffm in c("back","cent","fiad")) {
     for (ltt in 1:5)                    {
       numboxes<-numboxlist[ltt]
       Modelfunc(ltt,diffm,numboxes)
                                         }
                                           }  
ttext<-paste("Luminophore distribution at t = ",round(times[outtime],2),"year")
mtext(outer=TRUE,side=3,ttext,cex=2)