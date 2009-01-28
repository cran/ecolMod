###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 10.                               ##
## Dynamic programming                       ##
## Chapter 1.03. The patch selection model   ##
###############################################

#############################################################
## THE PATCH SELECTION MODEL                               ##
## Clark CW and M Mangel 1999.                             ##
## Dynamic State Variable Models in Ecology: methods and   ##
## applications. Oxford University Press.                  ##
## Chapter 1.                                              ##  
#############################################################


################################
# parameters settings          #
################################

x_crit  <- 0                 # critical mass to survive
x_max   <- 30                # maximal mass
x_rep   <- 4                 # critical mass for reproduction
x_class <- x_crit:x_max      # biomass classes
nmass   <- length(x_class)   # number of mass classes   

t_max   <- 20                # number of time steps 
times   <- 1:(t_max-1)
npatch  <- 3                 # number of patches

psurvive <- c(0.99,0.95,0.98)   # probability of surviving
pfood    <- c(0.2 ,0.5 ,0 )     # probability of feeding
cost     <- c(1 ,1 ,1 )         # cost of a patch
feedgain <- c(2 ,4 ,0 )         # gain of feeding
repr     <- c(0 ,0 ,4 )         # max reproduction

################################
# result matrices              #
################################

f         <- matrix(nrow=t_max,ncol=nmass ,0)    # optimal fitness values
bestpatch <- matrix(nrow=t_max-1,ncol=nmass-1,0) # best patch choice
V         <- vector(length=npatch)               # current fitness for a patch

# final fitness, at t_max
fend      <- 60
kx        <- 0.25*x_max

f[t_max,]  <- fend*(x_class-x_crit)/(x_class-x_crit+kx)

#####################################
# find fitness for mass x at time t #
#####################################

fitness <- function(x,t)

{
xx <- pmin(x ,x_max) 
xx <- pmax(xx,x_crit) 
fitness <- f[t,xx+1]
}


################################
# main optimisation loop       #
################################


for (t in rev(times))               # backward in time
{
 for (x in x_class[-1])             # for each biomass class, except x-crit
 {                                  
  dfit <- pmax(0,pmin(x-x_rep,repr)) # reproduction
  expectgain <- psurvive*( pfood   *fitness(x-cost+feedgain-dfit,t+1) + 
                          (1-pfood)*fitness(x-cost-dfit ,t+1) )
  V    <- dfit + expectgain  
  V[expectgain == 0] <- 0           # dead
  f[t,x+1]       <- max(V)          # optimal fitness
  bestpatch[t,x] <- which.max(V)    # best patch 
       
 }                                  # next biomass class x

}                                   # next time t

windows()
#par(mfrow=c(2,2))
cols <- c("black","darkblue","lightblue","white")
image(x=times,y=x_class[-1],z=bestpatch,ylab="weight",xlab="time",zlim=c(0,3),
      main="optimal patch",col=cols)
box()
legend("topleft",fill=cols,legend=c("dead","1","2","3"))
contour(x=1:t_max,y=x_class,z=f,add=TRUE)
legend("topright",legend="fitness",lty=1)
