###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 4.                                ##
## Parameterisation                          ##
## Chapter 4.4.1.                            ##
## Nonlinear fit of P-I data                 ##
###############################################

# the data
ll <- c(0.,1,10,20,40,80,120,160,300,480,700)
pp <- c(0.,1,3,4,6,8,10,11,10,9,8)

plot(ll,pp,xlab= expression("light, µEinst"~ m^{-2}~s^{-1}),
     ylab="production",pch=15,cex=1.5)

# Fitting
fit<-nls(pp ~pmax*2*(1+b)*(ll/iopt)/
                         ((ll/iopt)^2+2*b*ll/iopt+1),
         start=c(pmax=max(pp),b=0.005,iopt=ll[which.max(pp)]))

summary(fit)
pars <- as.list(coef(fit))

with(pars,
     curve(pmax*2*(1+b)*(x/iopt)/((x/iopt)^2+2*b*x/iopt+1),
     add=TRUE,lwd=2)
 )
title(expression (frac(2*pmax*(1+beta)*I/Iopt,
                 (I/Iopt)^2+2*beta*I/Iopt+1)),cex.main=0.8)

