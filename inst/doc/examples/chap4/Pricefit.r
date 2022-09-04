###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 4.                                ##
## Parameterisation                          ##
## Chapter 4.4.3.                            ##
## Pseudo-random search, a random-based      ##
## minimisation routine                      ##
###############################################

# the Price algorithm...
# included as a function in package ecolMod

pricefit <- function (
            par,                          # initial par estimates
            minpar=rep(-1e8,length(par)), # minimal parameter values
            maxpar=rep(1e8,length(par)),  # maximal parameter values
            func,                         # function to minimise
            npop=max(5*length(par),50),   # nr elements in population
            numiter=10000,                # number of iterations
            centroid = 3,                 # number of points in centroid 
            varleft  = 1e-8,              # relative variation upon stopping
            ...)

{

# Initialisation

 cost  <- function (par) func(par,...)
 npar  <- length(par)
 tiny  <- 1e-8
 varleft<-max(tiny,varleft)
 
 

 populationpar  <- matrix(nrow=npop,ncol=npar,byrow=TRUE,
             data= minpar+runif(npar*npop)*rep((maxpar-minpar),npop))
 colnames(populationpar)<-names(par)
 populationpar[1,]<-par

 populationcost <- apply(populationpar,FUN=cost,MARGIN=1)
 iworst         <- which.max(populationcost)
 worstcost      <- populationcost[iworst]


# Hybridisation phase   
 iter<-0
 while (iter<numiter & (max(populationcost)-min(populationcost))
                                  >(min(populationcost)*varleft))
 {
   iter<-iter+1
   
   selectpar <- sample(1:npop,size=centroid)  # for cross-fertilisation
   mirrorpar <- sample(1:npop,size=1)         # for mirroring  

   newpar    <- colMeans(populationpar[selectpar,])    # centroid
   newpar    <- 2*newpar - populationpar[mirrorpar,]   # mirroring
     
   newpar    <- pmin( pmax(newpar,minpar) ,maxpar)
                        
   newcost   <- cost(newpar)


    if (newcost < worstcost)  
     { 
       populationcost[iworst] <-newcost
       populationpar [iworst,]<-newpar  
       iworst     <- which.max(populationcost) # new worst member
       worstcost  <- populationcost[iworst]
     }
  } # end j loop


  ibest    <- which.min(populationcost)
  bestpar  <- populationpar[ibest,]
  bestcost <- populationcost[ibest]
return (list(par = bestpar, cost = bestcost, 
             poppar = populationpar, popcost=populationcost))

}

# The pseudo-data
amp    <- 6
period <- 5
phase  <- 0.5

x <- runif(20)*13 
y <- amp*sin(2*pi*x/period+phase) +rnorm(20,mean=0,sd=0.05)
plot(x,y,pch=16)


# Minimisation
# The model cost function
cost <- function(par) 
{
 with(as.list(par),{
    sum((amplitude*sin(2*pi*x/period+phase)-y)^2)
                   })
}

# Three methods tested
p1 <- optim(par=c(amplitude=1,period=1,phase=1), cost)
p2 <- optim(par=c(amplitude=1,period=1,phase=1), cost,method="SANN")
p3 <- pricefit(par=c(amplitude=1,period=1,phase=1),minpar=c(0,1e-8,0),
               maxpar=c(100,2*pi,100), func=cost,numiter=3000)

curve(p1$par[1]*sin(2*pi*x/p1$par[2]+p1$par[3]),lty=2,add=TRUE)
curve(p2$par[1]*sin(2*pi*x/p2$par[2]+p2$par[3]),lty=3,add=TRUE)
curve(p3$par[1]*sin(2*pi*x/p3$par[2]+p3$par[3]),lty=1,add=TRUE)

legend ("bottomright",lty=c(1,2,3),c("Price","Mathematical","Simulated annealing"))


