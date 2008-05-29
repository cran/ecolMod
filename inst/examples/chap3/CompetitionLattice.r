###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 3.                                ##
## Spatial components and transport          ##
## Chapter 3.6.4.                            ##
## Competition in a lattice grid -R-VERSION ##
###############################################

species <- c("Lolium","Agrostis","Holcus","Poa","Cynosurus")

replacement <- matrix(ncol=5,byrow=TRUE,data=c(
1.0,0.02,0.06,0.05,0.03,
0.23,1.0,0.09,0.32,0.37,
0.06,0.08,1.0,0.16,0.09,
0.44,0.06,0.06,1.0,0.11,
0.03,0.02,0.03,0.05,1.0 ) )


# competition function - THE SLOW VERSION
# See competition.f for the fast implementation in FORTRAN.

competition <- function(cells,nstep=100)

{
 nind  <- nrow(replacement)
 ncell <- nrow(cells)
 dens  <- matrix(nrow=nstep,ncol=nind)

 for (ss in 1:nstep)
  {
   dn  <- rbind(cells[ncell,]  ,cells[1:(ncell-1),])
   up  <- rbind(cells[2:ncell,],cells[1,]   )
   le  <- cbind(cells[,2:ncell],cells[,1]   )
   ri  <- cbind(cells[,ncell]  ,cells[,1:(ncell-1)])

   rnd <- matrix(nr=ncell,nc=ncell,runif(ncell*ncell))

   for (i in 1:ncell)
   {  for (j in 1:ncell)
     {  ii    <- cells[i,j]
        neigb <- c(up[i,j],ri[i,j],dn[i,j],le[i,j])
        p     <- replacement[neigb,ii]
        cump  <- c(cumsum(p/4),1)
        rep   <- min(which(cump>=rnd[i,j] ))
        if (rep<5) cells[i,j] <- neigb[rep]
      }
}
   for (i in 1:nind) dens[ss,i]<-sum(cells == i)
 }
return(list(cells=cells,density=dens))

}


##############################################
# application
##############################################

ini   <-c(4,5,1,3,2)
cells <- matrix(40,40,data=0)

cells[,1:8]  <-ini[1] ; cells[,9:16] <-ini[2] ;
cells[,17:24]<-ini[3] ; cells[,25:32]<-ini[4]
cells[,33:40]<-ini[5]

A100  <- competition(cells,nstep=100)
A200  <- competition(A100$cells,nstep=100)


# graphs

windows()
par(mfrow=c(2,2),mar=c(2,2,2,2))
col   <-c("grey","lightblue","blue","darkblue","black")

image(cells,col=col,zlim=c(1,5),axes=FALSE,main="initial")
box()
text(x=rep(0.1,5),y=seq(0.1,0.9,length.out=5),
     labels=species[ini],col=c("black",rep("white",4)),adj=0,font=2)

image(A100$cells,col=col, zlim=c(1,5),axes=FALSE,
      main="100 steps")
      box()
image(A200$cells,col=col, zlim=c(1,5),axes=FALSE,
      main="200 steps")
box()
matplot(rbind(A100$density,A200$density),type="l",lwd=2,lty=1,
        col=col,xlab="time",ylab="",axes=FALSE,frame.plot=TRUE)

