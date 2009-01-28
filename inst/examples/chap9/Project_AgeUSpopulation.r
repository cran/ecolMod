###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 9.                                ##
## Discrete time models                      ##
## Chapter 9.6.3 Project                     ##
## Equilibrium of US poppulation age class   ##
###############################################

##------------------------------------------------------------------------
# Keyfitz and Flieger 1971, from Caswell 2001
# US population of 1966, 5 year age classes, projection interval of 5 years
##------------------------------------------------------------------------

#------------------------#
# INITIALISE the model   #
#------------------------#


# Fecundity and Survival for each generation
NumClass       <- 10
Fecundity      <- c(0,0.00102,0.08515,0.30574,0.40002,0.28061,0.1526,0.0642,0.01483,0.00089)
Survival       <- c(0.9967,0.99837,0.9978,0.99672,0.99607,
                    0.99472,0.99240,0.98867,0.98274,NA)            # survival from i to i+1    
cbind(Fecundity,Survival)                                          # to check typos

# Population matrix M
DiffMatrix       <- matrix(data=0,nrow=NumClass,ncol=NumClass)     # declare it
DiffMatrix[1,]   <- Fecundity                                      # first row: fecundity
for (i in 1:(NumClass-1))  DiffMatrix[i+1,i] <- Survival[i]              

DiffMatrix                                                         # print the matrix to screen  

#------------------------#
# EIGENVALUE analysis:   #
#------------------------#

e <-eigen(DiffMatrix)
StableAge <- Re(e$vectors[,1])
StableAge <- StableAge/sum(StableAge)
StableAge                      # print
lambda    <- Re(e$values[1])
lambda                         # print 


# the reproductive value, uses the left eigen vectors...
lefte    <- eigen(t(DiffMatrix))
reprodVal<-Re(lefte$vectors[,1])
reprodVal<-reprodVal/reprodVal[1]  # a-dimensionalise
reprodVal                      # print

data.frame("class"=seq(2.5,47.5,by=5),"Stable"=StableAge,"Reprod"=reprodVal)

windows()
par(oma=c(0,0,3,0),mfrow=c(2,1))  
main <- paste ("Stable age distribution, rate of increase =",strtrim(lambda,5))
plot(StableAge,type="h",main=main,xlab="age class",ylab="-")
plot(reprodVal,type="h",main="Reproductive value",xlab="age class",ylab="-")
mtext(side=3,outer=TRUE,"Age-structured model",cex=1.5)
 
# Complex lambda
cbind(e$values,abs(e$values),atan(Im(e$values)/Re(e$values))/pi)

