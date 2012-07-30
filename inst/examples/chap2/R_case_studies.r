###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 2.                                ##
## Model formulation                         ##
## R case studies                            ##
###############################################

# Making sense out of mathematical formulations

r<-0.1
K<-10
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
curve (expr= r*x*(1-x/K),from=0,to=20,lwd=2,
       xlab="density, N",ylab="rate of change, dN/dt")
legend ("topright",c("K=10","r=0.1"))

# One formula, several parameter values
windows()
par(mar=c(5.1,4.1,4.1,2.1))
food <- seq(0,30,by=0.1)
ks   <- seq (1,10, by=2)

foodfun  <- outer(food,ks,function(x,y) x /(y +x))

matplot(x=food,foodfun,type="l",lty=1:10, col=1,
        xlab="food",ylab="-",
        main= expression (frac(food ,food+ks)))

legend ("bottomright", as.character(ks), title="ks=",
        col=1,lty=1:10)

