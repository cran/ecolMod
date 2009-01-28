###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 4.                                ##
## Parameterisation                          ##
## Chapter 4.4.2.                            ##
## Linear and Nonlinear fit of Pb210-data    ##
###############################################

# the data:
x    <- 0.5:9.5
y    <- c(3.9,1.7,1.1,0.5,0.3,0.2,0.1,0.05,0.03,0.02)

plot(y,x,pch=16,ylab="depth,cm", xlab="dpm/cm3",
     main="Pb210", ylim=c(10,0))

# the decay rate
lam  <- 0.031

# linear fit
LL   <- lm(log(y)~x)
C0   <- exp (coef(LL)[1])
Db   <- lam/(coef(LL)[2])^2

xx   <-seq(0,10,0.1)
lines(C0*exp(-sqrt(lam/Db)*xx),xx)

# nonlinear fit
fit<-nls(y ~C0*exp(-sqrt(lam/Db)*x),
         start=c(C0=5,Db=0.1))

C02  <- coef(fit)[1]
Db2  <- coef(fit)[2]

lines(C02*exp(-sqrt(lam/Db2)*xx),xx,lty=2)

legend("topleft",lty=c(1,2),
       c("linear fit", "nonlinear fit"))

par(new=TRUE,fig=c(0.5,1,0.1,0.6))

plot(y,x,pch=16,ylab="",xlab="",ylim=c(10,0),log="x")
lines(C0*exp(-sqrt(lam/Db)*xx),xx)
lines(C02*exp(-sqrt(lam/Db2)*xx),xx,lty=2)

