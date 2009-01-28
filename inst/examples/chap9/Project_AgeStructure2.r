###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 9.                                ##
## Discrete time models                      ##
## Chapter 9.6.2 Project                     ##
## Equilibrium dynamics of a simple          ##
## age-class model                           ##
###############################################

##------------------------------------------------------------------------
# population of long-lived individuals with 4 age classes
# example as in Gotelli, A Primer of Ecology. 
##------------------------------------------------------------------------

#------------------------#
# INITIALISE the model   #
#------------------------#


# Fecundity and Survival for each generation
AgeClasses       <- 4
ReprodRate       <- c(2,3,1,0)    
Survival         <- c(0.8,0.5,0.25,0.0)          # survival from i to i+1

# Population matrix M
Leslie       <- matrix(data=0,nrow=AgeClasses,ncol=AgeClasses) # declare it
Leslie[1,]   <- ReprodRate * Survival[]                     # first row :fecundity
for (i in 1:(AgeClasses-1))  Leslie[i+1,i] <- Survival[i]              

Leslie                   # print the matrix to screen      

#------------------------#
# RUN the model          #
#------------------------#

numsteps    <- 50
Population  <- rep(0.1,times=AgeClasses) # initial condition
out         <- Population  

for (i in 1:numsteps)
{
 Population <- Leslie %*%Population        # %*% is matrix multiplication
 out        <- rbind (out, t(Population))  # row bind: add output row to matrix 'out'
}

#------------------------#
# PLOTTING model output: #
#------------------------#

# First a function to plot histogram-like figures is defined

Ageplot <- function (Age,main=NULL,xlab=NULL,ylab=NULL,axes=TRUE,fill=NULL)
 {
N <- length(Age)
y <- vector(length=2*N)
x <- vector(length=2*N)
for (i in 1:N) y[(2*i-1):(2*i)]<-Age[i]
y <- c(0,y,0)
for (i in 1:N) x[(2*i-1):(2*i)]<-i
x <- c(0,0,x)

plot(x,y,type="l",main= main,xlab=xlab,ylab=ylab,xlim=c(0,N),axes=axes)
if (!is.null(fill)) polygon(x,y,col=fill)
lines(1:N,Age,type="h")
if (!axes)box(col="grey")

}

# 1. the age classes 

windows()                            
par(mfrow=c(4,4), oma=c(0,0,3,0),mar=c(1,1,3,0))   # set number of plots (mfrow) and outer and inner margin 
for (i in 1:16)Ageplot(out[i,],main=paste("step",i),xlab="",ylab="",axes=FALSE,fill="grey")

mtext(side=3,outer=TRUE,"Age-structured model",cex=1.5)

# 2. The changes in time 

windows(width=8.5,height=11)                       ## A4
par(mfrow=c(3,2), oma=c(0,0,3,0))                  # set number of plots (mfrow) and margin size (oma)

# AgeClass densities versus time
plot (0:numsteps,out[,1],type="l",main="ageclass densities" ,xlab="time, days",ylab="-",log="y")
for (i in 2:AgeClasses) lines(0:numsteps,out[,i],lty=i)
legend("topleft",legend=1:4,lty=1:4)

# total density versus time
Density  <- vector(length=numsteps+1)              
for (i in 1: numsteps+1) Density  [i] <- sum(out[i,])
plot (0:numsteps,Density,type="l",main="total density"    ,xlab="time, days",ylab="-")

# mean age versus time
MeanAge   <- vector(length=numsteps+1)             
for (i in 1: numsteps+1) MeanAge[i] <- out[i,]%*%seq(1:AgeClasses)/sum(out[i,])
plot (0:numsteps,MeanAge,type="l",main="mean age"    ,xlab="time, days",ylab="-")

# Rate of increase: ci+1/ci
Increase  <- Density[2:(numsteps+1)]/Density[1:numsteps] 
plot (1:numsteps,Increase,type="l",main="finite rate of increase"    ,xlab="time, days",ylab="-")

# Final age distribution 
FinalAge <- out[numsteps+1,]                        
FinalAge <- FinalAge /sum(FinalAge)
Ageplot(FinalAge ,main="Final age distribution",xlab="age",ylab="-")

mtext(side=3,outer=TRUE,"Age-structured model",cex=1.5)


#------------------------#
# EIGENVALUE analysis:   #
#------------------------#

# the stable age distribution, uses the right eigen vectors...
e <-eigen(Leslie)
StableAge <- Re(e$vectors[,1])
StableAge <- StableAge/sum(StableAge)
StableAge                      # print
FinalAge                       # and compare..

# the rate of increase, the largest eigenvalue
lambda    <- Re(e$values[1])
rincrease <- log(lambda)
c(lambda,rincrease)            # print 
Increase[numsteps]             # and compare

# the reproductive value, uses the left eigen vectors...
lefte    <- eigen(t(Leslie))
reprodVal<-Re(lefte$vectors[,1])
reprodVal<-reprodVal/reprodVal[1]  # a-dimensionalise

windows()
par(oma=c(0,0,3,0),mfrow=c(2,1))  
main <- paste ("Stable age distribution, rate of increase =",strtrim(lambda,5))
Ageplot(StableAge,main=main,xlab="age",ylab="-",fill="grey")
Ageplot(reprodVal,main="Reproductive value",xlab="age",ylab="-",fill="grey")
mtext(side=3,outer=TRUE,"Age-structured model",cex=1.5)
 

