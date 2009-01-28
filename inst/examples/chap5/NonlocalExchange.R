###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 5.                                ##
## Model solution - analytical methods       ##
## Chapter 5.4.4                             ##
## Nonlocal exchang sediment model           ##
###############################################

#-------------------------------------------------------------------*
# Nonlocal exchange model                                           *
# concentration of a substance subjected to decay (rate coeff k),   *
# to diffusive mixing (diff coeff Db) and to nonlocal exchange, i.e.*
# part of the flux (p) directly injected in the sediment at depth L *
# EQUATIONS:                                                        *
# [0-L]: sedimentation, mixing, consumption                         *
# dC/dt=0 = - w dC/dx + Db (d2C/dx2) - k C                          *
# [L-Inf]: sedimentation, mixing, consumption                       *
# dC/dt=0 = - w dC/dx + Db (d2C/dx2) - k C                          *
# Solving this system gives an equation with 3 integration constants*
# A1, A2, A3, which are defined by the boundary conditions          *

# BOUNDARIES:                                                       *
# x=0, surface                                                      *
#    deposition Flux : w C|x=0 -Db (dC/dx)|x=0 = flux               *
# injection depth L                                                 *
#    Continuity of concentration: conc1(L) = conc2(L)               *
#    Continuity of flux,
# injection flux is added to advective/diffusive flux of layer 1 *
# This gives a set of 3 linear equations in the 3 unknown           *
# integration constants, which is solved by matrix inversion        *
# After solving this parameter, the vertical profile is calculated  *
#-------------------------------------------------------------------*


nonlocal <- function( 
        k=0.03108,v=0.001,Db=0.1,depo=0.1,L=5,injectflux=0.3,sed
                     )
{
 # Power and exponents
 a  <- (v/Db - sqrt ( (v/Db)^2 + 4*k / Db))/2.
 b  <- (v/Db + sqrt ( (v/Db)^2 + 4*k / Db))/2.
 expaL <- exp(a*L)
 expbL <- exp(b*L)
# coefficients multiplying with integration cts 
 A <- matrix(nrow=3,ncol=3,byrow=TRUE,  data=   c(
#  col1: coeffA1  , col2: coeff A2, col3: coeff A3
   v-Db*a         ,v-Db*b         ,0             , 
   expaL          ,expbL          ,-expaL        ,   
(v-Db*a)*expaL  ,(v-Db*b)*expbL ,(Db*a-v)*expaL )
           )
# right hand side
#      Surface flux, continuity conc, continuity flux
 B  <- c(depo     ,  0            ,- injectflux)
# solve
 X  <- solve(A,B)

 s1 <- which (sed<L)
 s2 <- which (sed>=L)

 conc <- vector(length=length(sed))
 conc[s1] <- X[1]* exp(a*sed[s1])+X[2]*exp(b*sed[s1])
 conc[s2] <- X[3]* exp(a*sed[s2])
 return(conc)

}   # end nonlocal


#-------------------------------------------------------------------*
# Model application 1: Pb210
#-------------------------------------------------------------------*

depth <- seq(0,15,by=0.01)              # cm          depth of sediment slices

# First set of parameters
k            <- 0.03108                 # /yr         first-order decay rate
Db           <- 0.1                     # cm2/yr      bioturbation coefficient
w            <- 0.001                   # cm/yr       advection rate
depo         <- 0.1                     # dpm/cm2/yr  Pb210 deposition flux
injecflux    <- 0.3                     # dpm/cm2/yr  Pb210 injection flux
injecdepth   <- 5                       # cm          Pb210 injection depth

Pb <-nonlocal(k,w,Db,depo,injecdepth,injecflux,sed=depth)
plot(Pb,depth,ylim=c(15,0),xlab="Pb, dpm/cm3",ylab="cm",type="l",lwd=2,main="75% injected")

