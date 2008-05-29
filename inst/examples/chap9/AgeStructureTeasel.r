###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 9.                                ##
## Discrete time models                      ##
## Chapter 9.6                               ##
## Population dynamics of Teasel             ##
###############################################

##------------------------------------------------------------------------
# Teasel example from Caswell 1989
# Teasel, a perennial weed is modeled in six stages
##------------------------------------------------------------------------

#------------------------#
# INITIALISE the model   #
#------------------------#

# Define the stages
Stagenames <-c("dormant seeds 1yr","dormant seeds 2yr",
"rosettes<2.5cm","rosettes2.5-18.9cm","rosettes>19cm",
"flowering plants")

NStages   <- length(Stagenames)
# Population matrix M

M <- matrix(nrow=NStages,ncol=NStages,byrow=TRUE,data = c(
    0,      0,      0,      0,      0,      322.38,
    0.966,  0,      0,      0,      0,      0     ,
    0.013,  0.01,   0.125,  0,      0,      3.448 ,
    0.007,  0,      0.125,  0.238,  0,      30.170,
    0.008,  0,      0,      0.245,  0.167,  0.862 ,
    0,      0,      0,      0.023,  0.75,   0      )  )

rownames(M) <- Stagenames
colnames(M) <- 1:NStages

M   

#------------------------#
# Stable stage + lambda  #
#------------------------#

numsteps      <- 50
Population    <- rep(0,times=NStages)
Population[1] <- 100
out           <- matrix(nrow=numsteps+1,ncol=NStages,data=0) 
out[1,]       <- Population

for (i in 1:numsteps)
{
 Population <- M %*%Population    
 out[i+1,]  <- t(Population)  
}

#------------------------#
# PLOTTING  
#------------------------#

windows(5,5)
p   <- out/rowSums(out)
matplot(p,type="l",lty=1:6,main="stage distribution",col="black")
legend("topright", legend=Stagenames,lty=1:6)

p[nrow(p),]
sum(out[50,])/sum(out[49,])


e <-eigen(M)
Stablestage <- Re(e$vectors[,1])
Stablestage <- Stablestage/sum(Stablestage)
                       
lambda    <- Re(e$values[1])

#------------------------#
# The reproductive value #
#------------------------#

rPop     <- diag(nrow=NStages,ncol=NStages)

for (i in 1:100)
  for (j in 1:NStages) rPop [j,] <- M %*% rPop [j,]    

rr<-rowSums(rPop)
rr<-rr/rr[1] 

lefte     <- eigen(t(M))
reprodVal <-Re(lefte$vectors[,1])
reprodVal <-reprodVal/reprodVal[1]  

rr
reprodVal


lambda
data.frame("stage"=Stagenames,"Stable"=Stablestage,"Reprod"=reprodVal)
