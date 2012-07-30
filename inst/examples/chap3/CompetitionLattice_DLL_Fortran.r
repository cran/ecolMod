###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 3.                                ##
## Spatial components and transport          ##
## Chapter 3.6.4.                            ##
## Competition in a lattice grid             ##
## Function in DLL programmed in Fortran     ##
###############################################

# First compile the fortran source:
#wd<-getwd()
#setwd(" the working directory here")
#system("R CMD SHLIB competition.f")

# load the DLL
dll <- "competition.dll"
dyn.load(dll)


# competition function in FORTRAN
competition <- function(cells,nstep=100)

{
 nspec   <- nrow(replacement)
 ncell   <- nrow(cells)
 cells   <- matrix(nrow=ncell,ncol=ncell,data=as.integer(cells))
 sumdens <- matrix(nrow=nstep,ncol=nspec,as.integer(0))
 seed    <- runif(1)*-10
 res <- .Fortran("lattice",nspec=nspec,ncell=ncell,nstep=as.integer(nstep),
                 cells=cells,replacement=as.double(replacement),
                 sumdens=sumdens,seed=as.integer(seed))
 
 return(list  (cells=matrix(nr=ncell,nc=ncell,res$cells),
               density=res$sumdens))

}


##############################################
# application
##############################################
species <- c("Lolium","Agrostis","Holcus","Poa","Cynosurus")

replacement <- matrix(ncol=5,byrow=TRUE,data=c(
1.0,0.02,0.06,0.05,0.03,
0.23,1.0,0.09,0.32,0.37,
0.06,0.08,1.0,0.16,0.09,
0.44,0.06,0.06,1.0,0.11,
0.03,0.02,0.03,0.05,1.0 ) )


ini   <-c(4,5,1,3,2)
cells <- matrix(40,40,data=0)

cells[,1:8]  <-ini[1] ; cells[,9:16] <-ini[2] ;
cells[,17:24]<-ini[3] ; cells[,25:32]<-ini[4]
cells[,33:40]<-ini[5]

nstep=100
A100  <- competition(cells,nstep=100)
A200  <- competition(A100$cells,nstep=100)


# graphs

par(mfrow=c(2,2),mar=c(2,2,2,2))
col   <-c("grey","lightblue","blue","darkblue","black")

image(cells,col=col,zlim=c(1,5),axes=FALSE,main="initial")
text(x=rep(0.1,5),y=seq(0.1,0.9,length.out=5),
     labels=species[ini],col="white",adj=0,font=2)
image(A100$cells,col=col, zlim=c(1,5),axes=FALSE,
      main="100 steps")
image(A200$cells,col=col, zlim=c(1,5),axes=FALSE,
      main="200 steps")
matplot(rbind(A100$density,A200$density),type="l",lwd=2,lty=1,
        col=col,xlab="time",ylab="",axes=FALSE,frame.plot=TRUE)

dyn.unload(dll)
