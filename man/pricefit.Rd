\name{pricefit}
\alias{pricefit}
\title{Pseudo-random search algorithm of Price (1997)}
\description{Pseudo-random search algorithm of Price (1997). Usesd in the book as an example of 
a random-based fitting technique, and as an example of how to create a function in R.}
\usage{pricefit(par, minpar=rep(-1e8,length(par)),
maxpar=rep(1e8,length(par)), func, npop=max(5*length(par),50),
numiter=10000, centroid = 3, varleft  = 1e-8,...)}
\arguments{
  \item{par }{initial values of the parameters to be optimised}
  \item{minpar }{minimal values of the parameters to be optimised}
  \item{maxpar }{maximal values of the parameters to be optimised}
  \item{func }{function to be minimisedA function to be minimized, its first argument should bw the vector of parameters
    over which minimization is to take place. It should return a scalar result, the model cost, e.g the sum of squared residuals.}
  \item{npop }{number of elements in population}
  \item{numiter }{number of iterations to be performed. Defaults to 10000. There is no other stopping criterion. }
  \item{centroid }{number of elements from which to estimate new parameter vector}
  \item{varleft }{relative variation remaining; if below this value algorithm stops}
  \item{... }{arguments passed to funtion func}
}
\value{
  a list containing:
  \item{par }{the optimised parameter values}
  \item{cost }{the model cost, or function evaluation associated to the optimised parameter values, i.e. the minimal cost}
  \item{poppar }{all parameter vectors remaining in the population, matrix of dimension (npop,length(par))}
  \item{popcost }{model costs associated with all population parameter vectors, vector of length npop}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{pricefit  # will display the code 

amp    <- 6
period <- 5
phase  <- 0.5

x <- runif(20)*13 
y <- amp*sin(2*pi*x/period+phase) +rnorm(20,mean=0,sd=0.05)
plot(x,y,pch=16)


cost <- function(par) sum((par[1]*sin(2*pi*x/par[2]+par[3])-y)^2)

p1 <- optim(par=c(amplitude=1,phase=1,period=1), cost)
p2 <- optim(par=c(amplitude=1,phase=1,period=1), cost,method="SANN")
p3 <- pricefit(par=c(amplitude=1,phase=1,period=1),minpar=c(0,1e-8,0),
               maxpar=c(100,2*pi,100), func=cost,numiter=3000)

curve(p1$par[1]*sin(2*pi*x/p1$par[2]+p1$par[3]),lty=2,add=TRUE)
curve(p2$par[1]*sin(2*pi*x/p2$par[2]+p2$par[3]),lty=3,add=TRUE)
curve(p3$par[1]*sin(2*pi*x/p3$par[2]+p3$par[3]),lty=1,add=TRUE)
legend ("bottomright",lty=c(1,2,3),c("Price","Mathematical","Simulated annealing"))}
\details{see the book of Soetaert and Herman for a description of the algorithm AND for a line to line explanation of the function code.}
\seealso{\code{\link{optim}} the R default.}
\keyword{optimize}

