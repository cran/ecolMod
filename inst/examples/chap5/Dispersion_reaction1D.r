###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 5.                                ##
## Model solution - analytical methods       ##
## Chapter 5.4.1. and 5.4.2                  ##
## Transient dispersion reaction in one      ##
## dimension and on a 2-dimensional surface  ##
###############################################

require(shape)  # colors, drapecol...

##########################################
# The diffusion-reaction equation in 1-D #
# cartesian and cylindrical coordinates  #
##########################################

Ds  <- 1    # diffusion coefficient
ini <- 1    # initial condition
k   <- 0.05 # growth rate

#=============================
# 1-D cartesian coordinates  
#=============================
windows(5,5)
grow1D <- outer(xx<-seq(-5,5,length=50),tt<-seq(0.1,5,by=0.025),
                FUN = function (x,tt) ini/(2*sqrt(pi*Ds*tt))*exp(k*tt-x^2/(4*Ds*tt)))

persp(xx,tt,z=grow1D,theta=150,box=TRUE,axes=TRUE,col=drapecol(grow1D,femmecol(100)),
xlab="space",ylab="time",zlab="Conc",border=NA,main="1-D diffusion-reaction")  


#=============================
# 1-D cylindrical coordinates  
#=============================
                         
growplane <- outer(rr<-seq(-5,5,length=50),tt<-seq(0.1,5,by=0.025),
              FUN = function (rr,tt) ini/(4*pi*Ds*tt)*exp(k*tt-rr^2/(4*Ds*tt)))

plotplane <- function(time,rmax=5,...) 
 { 
  val <-outer(xx<-seq(-rmax,rmax,length=50),yy<-xx,
              FUN = function (x,y) 
             {
              r2<-x*x+y*y;
              ini/(4*pi*Ds*time)*exp(k*time-r2/(4*Ds*time))
              } 
             ) 
  persp(xx,yy,z=val,theta=150,box=TRUE,axes=TRUE,col=drapecol(val,femmecol(100)),
        zlab="conc",border=NA,...)#,border=NA)          
 }

windows()
par(mfrow=c(2,2),mar=c(3,3,3,3))
plotplane(0.1, main= "0.1 day")
plotplane(1  , main= "1 day")
plotplane(2  , main= "2 days") 
plotplane(5  , main= "5 days")  



