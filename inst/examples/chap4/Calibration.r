###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 4.                                ##
## Parameterisation                          ##
## Chapter 4.4.4.                            ##
## Calibration of a simple model             ##
###############################################

# this takes a long time...
require(ecolMod)   # knows pricefit

# Function to calculate model cost  

costf <- function(params)
  {with(as.list(params),
   {
      Carbon    <- meanDepo*mult/k
      outtimes  <- as.vector(oxcon$time)
      outmin    <- ode(Carbon,outtimes,minmod,params)
      costt     <- sum((outmin[,3]-oxcon$cons)^2)
      return (costt)
   })
  }

# Function to calculate rate of change of state variable

minmod <- function(t,Carbon,parameters)
 {  with (as.list(c(Carbon,parameters)),
     {
        minrate  <- k*Carbon
        Depo     <- approx(Flux[,1],Flux[,2], xout=t)$y
        dCarbon  <- mult*Depo - minrate
    list(dCarbon,minrate)
     })
 }
 
 
 
require(deSolve)
 
# Define problem and data
 
Flux <- matrix(ncol=2,byrow=TRUE,data=c(
  1, 0.654, 11, 0.167, 21, 0.060, 41, 0.070,
 73, 0.277, 83, 0.186, 93, 0.140,103, 0.255,
113, 0.231,123, 0.309,133, 1.127,143, 1.923,
153, 1.091,163, 1.001,173, 1.691,183, 1.404,
194, 1.226,204, 0.767,214, 0.893,224, 0.737,
234, 0.772,244, 0.726,254, 0.624,264, 0.439,
274, 0.168,284, 0.280,294, 0.202,304, 0.193,
315, 0.286,325, 0.599,335, 1.889,345, 0.996,
355, 0.681,365, 1.135))

meanDepo   <- mean(approx(Flux[,1],Flux[,2], xout=seq(1,365,by=1))$y)

oxcon<-as.data.frame(matrix(ncol=2,byrow=TRUE,data=c(
 68, 0.387, 69, 0.447, 71, 0.473, 72, 0.515,
189, 1.210,190, 1.056,192, 0.953,193, 1.133,
220, 1.259,221, 1.291,222, 1.204,230, 1.272,
231, 1.168,232, 1.168,311, 0.963,312, 1.075,
313, 1.023)))

names(oxcon)<-c("time","cons")


multser   <- seq(1,1.5,by=.05)
numms     <- length(multser)
kseries   <- seq(0.001,0.05,by=0.002)
numks     <- length(kseries)

outcost <- matrix(nrow=numms,ncol=numks)

    for (m in 1:numms)
    {
    for (i in 1:numks)
     {
      pars         <- c(k=kseries[i],mult=multser[m])
      outcost[m,i] <- costf(pars)
      }
     }

minpos<-which(outcost==min(outcost),arr.ind=TRUE)
multm<-multser[minpos[1]]
ki<-kseries[minpos[2]]

# 3 runs...
optpar <- pricefit(par=c(k=ki,mult=multm),minpar=c(0.001,1),
          maxpar=c(0.05,1.5),func=costf,npop=50,numiter=500,
          centroid=3,varleft=1e-8)

optpar20 <- pricefit(par=optpar$par,minpar=c(0.001,1),
            maxpar=c(0.05,1.5),func=costf,npop=50,numiter=500,
            centroid=3,varleft=0.2)
optpar25 <- pricefit(par=optpar$par,minpar=c(0.001,1),
            maxpar=c(0.05,1.5),func=costf,npop=50,numiter=500,
            centroid=3,varleft=0.025)


windows()
outtimes      <- seq(1,365,by=1)
Carbon        <- meanDepo*optpar$par[2]/optpar$par[1]
names(Carbon) <-"Carbon"
out           <- as.data.frame(ode(Carbon,outtimes,minmod,optpar$par))
names(out)    <- c("time","Carbon","minrate")

par (oma=c(0,0,0,2))
plot(Flux,type="l",xlab="daynr",ylab="mmol/m2/d",
     main="Sediment-detritus model",lwd=2)
lines(out$time,out$minrate,lwd=2,col="darkgrey")
points(oxcon$time,oxcon$cons,pch=25,col="black",bg="darkgray",cex=2)
par(new=TRUE)
plot(out$time,out$Carbon,axes=FALSE,xlab="",ylab="",
       type="l",lty=2)
axis(4)
mtext(side=4,"mmolC/m2",outer=TRUE)
legend("topleft",col=c("black","darkgrey","black"),
leg=c("C flux","C mineralisation","C concentration"),
      lwd=c(2,2,1),lty=c(1,1,2))
        

windows()
pgr<-gray.colors(n=25, start=0.95, end=0.0)
filled.contour(x=multser,y=kseries,z=outcost,
             ylab="k (/day)",xlab="multiplication factor (-)",
             main="Model cost landscape",col=pgr,nlevels=25,
             plot.axes={ 
             axis(1); axis(2);
             points(optpar20$poppar[,2],optpar20$poppar[,1],pch="o",cex=.5);
             points(optpar25$poppar[,2],optpar25$poppar[,1],pch="+",cex=1);
             points(optpar$par[2],optpar$par[1],pch=16,cex=2) 
                       }
               )

