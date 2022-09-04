###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 11.                               ##
## Testing and validating the model          ##
## Chapters 11.6.3                           ##
## Univariate and bivariate sensitvity       ##
## The Bacteria+glucose model                ##
###############################################

require(deSolve)

model <- function(t,state,pars)
{
with (as.list(c(state,pars)), {

dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact

return(list(c(dBact,dSub)))
                              })
}

pars <- list(Bini=0.1,Sini=100,gmax =0.5,eff = 0.5,
              ks =0.5, rB =0.01, dB =0.01)

tout    <- seq(0,50,by=0.5)
state   <- c(Bact=pars$Bini,Sub =pars$Sini)
out     <- as.data.frame(ode(state,tout,model,pars))

plot(out$time,out$Bact,ylim=range(c(out$Bact,out$Sub)),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
lines(out$time,out$Sub,lty=2,lwd=2)
lines(out$time,out$Sub+out$Bact)

legend("topright",c("Bacteria","Glucose","TOC"),
       lty=c(1,2,1),lwd=c(2,2,1))



windows()
Reference <- out
pp        <- unlist(pars)
tiny      <- 0.1
par(mfrow=c(3,3))
for (i in 1:length(pars))
{
 pars[i] <- pp[i]*(1+tiny)
 state   <- c(Bact=pars$Bini,Sub =pars$Sini)
 out     <- as.data.frame(ode(state,tout,model,pars))

 plot(out$time,out$Bact,xlab="hour",ylab="molC/m3",
      type="l",lwd=2,main=names(pars)[i])
 lines(Reference$time,Reference$Bact,lty=2)
 pars[i] <- pp[i]
}
plot(0,axes=FALSE,xlab="",ylab="",type="n")
legend("center",c("perturbed","reference"),lwd=c(2,1),lty=c(1,2))


###############################
# Fig. 11.7 Sensitivity functions
###############################

windows()
state <- c(Bact=pars$Bini, Sub =pars$Sini)
yRef  <- as.data.frame(ode(state,tout,model,pars))$Bact
pp    <- unlist(pars)
nout  <- length(yRef)
npar    <- length(pars)

# perturbation factor
tiny    <- 1e-8
dp      <- pp*tiny

Sens    <- matrix(nrow=nout,ncol=npar,NA)

for (i in 1:npar)
{
  dval    <- pp[i]+dp[i]
  pars[i] <- dval
  state   <- c(Bact=pars$Bini,Sub =pars$Sini)
  yPert   <- as.data.frame(ode(state,tout,model,pars))$Bact
  Sens[,i]<- (yPert-yRef)/tiny
  pars[i] <- pp[i]
}
colnames(Sens) <- names(pars)
rownames(Sens) <- tout
format(as.data.frame(Sens[1:5,]),digits=2)

par(mfrow=c(1,1))
matplot(tout,Sens,type="l",lty=1 :10,col=1,lwd=2,main="Sensitivity functions")
legend("topright",names(pars),lty=1:10,col=1,lwd=2)


mabs <- colMeans(abs(Sens))
msqr <- sqrt(colSums(Sens*Sens)/nout)
format(data.frame(msqr,mabs),digits=2)

###############################
# Fig. 11.8 bivariate analysis 
###############################

windows()
#par(oma=c(3,3,3,3))
panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))

pairs(Sens,upper.panel=panel.cor)
mtext(outer=TRUE,side=3,line=1.5,
      "Sensitivity functions",cex=1.5)
