###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 5.                                ##
## Model solution - analytical methods       ##
## Chapter 5.5.3. Project                    ##
## Carbon dynamics in sediment               ##
###############################################

#================================
# The analytical model
#================================

Cprofile <- function(k,        # /day        first-order decay rate
                     Db,       # cm2/d       mixing (bioturbation) coefficient
                     w,        # cm/d        advection
                     flux,     # mmol/cm2/d  carbon deposition flux
                     mix,      # cm          mixing depth 
                     depth)    # cm          vector with depth of layers
{

# Generates a carbon-sediment depth profile #

  # three exponents: a1,b1 for first layer {C=A1*e^(a1*x)+B1*e^(b1*x)}
  a1    <- (w-sqrt(w*w+4*Db*k))/2/Db
  b1    <- (w+sqrt(w*w+4*Db*k))/2/Db
  #                : a2, for second layer {C= A2*e^(a2*x)}
  a2    <- -k/w

  # three equations in 3 unknown integration constants (A1,B1,A2)
  # equations are written as: A*X = B
  A   <- matrix (nrow=3,ncol=3)
  B   <- vector (length=3)

  # collect terms for the unknowns
  # eq 1: flux equation: flux=-A1*a1*Db - B1*b1*Db + w*A1 + w*B1
  A[1,] <- c(-a1*Db + w , - b1*Db + w, 0)
  B[1]  <- flux

  # eq 2: concentration continuity at depth mix:
  # A1*e^(a1*mix) + B1*e^(b1*mix) - A2* e^(a2*mix)=0
  A[2,] <- c(exp(a1*mix) ,exp(b1*mix), -exp(a2*mix))
  B[2]  <- 0

  # eq 3: flux continuity at depth mix:
  #-A1*a1*Db*e^(a1*mix) - B1*b1*Db*e^(b1*mix) +w*A1*e^(a1*mix) + w*B1*e^(b1*mix)
  #                                      - w*A2*e^(a2*mix)=0

  A[3,] <- c((-a1*Db+w)*exp(a1*mix) ,(-b1*Db+w)*exp(b1*mix), -w*exp(a2*mix))
  B[3]  <- 0

  # Try to Solve linear equation for A1, B1, A2 
  X      <- try(solve(A,B),silent=TRUE)

  # Calculate Carbon in two sections of sediment (bioturbated/not bioturbated)
  Carbon <- vector (length=length(depth))  
  
  if (length(X) == 3)   {  # Solution was found
  c1 <- which(depth<mix) 
  c2 <- which(depth>=mix)

  Carbon[c1] <- X[1]*exp(a1*depth[c1])+X[2]*exp(b1*depth[c1])
  Carbon[c2] <- X[3]*exp(a2*depth[c2])
                         }
  # If not succesfull: A is singular, no C escapes the mixed layer
  #                    the solution is simple
  if (length(X)!= 3) Carbon <-flux/(w-a1*Db)*exp(a1*depth)

  return (Carbon)

 } # end function carbon profile

#================================
# The model application
#================================

windows(width=8.5,height=11)         ## A4
par(oma = c(0,0,3,0), mfrow=c(2,2))  ## outer margin, multiple figures 2 rows, 2 cols

depth <- seq(0,10,by=0.01)           # cm          depth of slices into sediment

# First set of parameters
k     <- 0.01                        # /day        first-order consumption rate
Db    <- 1/365                       # cm2/day     diffusion coefficient
w     <- 1/365                       # cm/day      advection rate
flux  <- 0.1                         # mmol/cm2/d  carbon deposition flux
mix   <- 3                           # cm          mixed layer depth    

Car <-Cprofile (k,Db,w,flux,mix,depth)
plot(Car,depth,ylim=c(10,0),xlab="C,mmol/cm3",ylab="cm",type="l",lwd=2,main="k=0.01/d,Db=1cm2/yr")

# Sensitivity analyses
krange <- c(0.002,0.004,0.006,0.008,0.01)          # k-values
Cmat   <- NULL                                     # matrix with profiles
for (k_ in krange) Cmat<-cbind(Cmat,Cprofile(k_,Db,w,flux,mix,depth))
matplot(y=depth,x=Cmat,ylim=c(10,0),xlab="C,mmol/cm3",ylab="cm",type="l",lty=1,lwd=2)
legend ("bottomright",as.character(krange),lty=1,lwd=2,col=1:10,title= "k, /day")


Dbrange <- c(0.001,0.01,0.1,1,10)/365                # Db-values, cm2/day
Cmat   <- NULL                                       # matrix with profiles
for (Db_ in Dbrange) Cmat<-cbind(Cmat,Cprofile(k,Db_,w,flux,mix,depth))

matplot(y=depth,x=Cmat,ylim=c(10,0),xlab="C,mmol/cm3",ylab="cm",type="l",lty=1,lwd=2)
legend ("bottomright",as.character(Dbrange*365),lty=1,lwd=2,col=1:10,title= "Db, cm2/yr")


wrange <- c(0.001,0.01,0.1,1,10)/365                 # w-values, cm/day
Cmat   <- NULL                                       # matrix with profiles
for (w_ in wrange) Cmat<-cbind(Cmat,Cprofile(k,Db,w_,flux,mix,depth))

matplot(y=depth,x=Cmat,ylim=c(10,0),xlab="C,mmol/cm3",ylab="cm",type="l",lty=1,lwd=2)
legend ("bottomright",as.character(wrange*365),lty=1,lwd=2,col=1:10,title= "w, cm/yr")

mtext(outer=TRUE,side=3,"Analytical model - Org C in sediment",cex=1.5)

