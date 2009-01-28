###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 6.                                ##
## Model solution - numerical methods        ##
## Chapter 6.6.3                             ##
## 0-D estuarine zooplankton model           ##
###############################################

#-----------------------#
# the forcing function: #
#-----------------------#
fZooTime = c(0, 30,60,90,120,150,180,210,240,270,300,340,367)
fZooConc = c(20,25,30,70,150,110, 30, 60, 50, 30, 10, 20, 20)

# The model, integrated with EULER
euler <-function(start,end,delt,g,k)
 {
  times     <- seq(start,end,delt)
  nt        <- length(times)

  ZOOsea    <- approx(fZooTime,fZooConc, xout=times)$y

  ZOO       <- 5
  out       <- matrix(ncol=3,nrow=nt)

  for (i in 1:(nt-1))
   {
     decay   <- g*ZOO
     input   <- k*(ZOOsea[i] - ZOO)

     out[i,] <- c(times[i],ZOO, input)

     dZOO    <- input - decay

     ZOO     <- ZOO + dZOO *delt
   }

  colnames(out) <- c("time","ZOO","input")
  return(as.data.frame(out))

}

out <-euler(0,365,0.1,k=0.015,g=0.05)

par (oma=c(0,0,0,2))

plot(fZooTime,fZooConc,type="b",xlab="daynr",
     ylab="gC/m3",pch=15,lwd=2, 
     main="Zooplankton model")
lines(out$time,out$ZOO,lwd=2,col="darkgrey")

legend("topright",c("Marine zooplankton","Estuarine zooplankton"),
       pch=c(15,NA),lwd=2,col=c("black","grey"))
