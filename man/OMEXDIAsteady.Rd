\name{OMEXDIAsteady}
\alias{OMEXDIAsteady}
\title{steady-state application of the OMEXDIA diagenetic model - fortran DLL}
\description{A 1-D model of Carbon, nitrogen and oxygen diagenesis in a marine sediment.
\cr The model describes six state variables, in \bold{100} layers:
  \itemize{
\item 2 fractions of organic carbon (FDET,SDET): fast and slow decaying, solid substance
\item Oxygen (O2), dissolved substance
\item Nitrate (NO3), dissolved substance
\item Ammonium (NH3), dissolved substance
\item Oxygen demand unit (ODU), dissolved substance, as lump sum of all reduced substances other than ammonium
}
See Soetaert et al., 1996 for further details of the model.

This is a simplified version of the OMEXDIA model, added just to create a figure in the book.
\cr A more complete version will be published in a separate R-package that will deal with reactive transport modelling in R.
\cr The name of this package is not yet decided upon.}
\usage{OMEXDIAsteady()}
\value{
  a list containing:
  \item{steady }{The steady-state condition of the state variables, a vector containing steady-state concentrations of
FDET(0-100), SDET(101-200), O2 (201-300), NO3 (301-400), NH3 (401-500) and ODU (501-600)}
  \item{precis }{the precision of the steady-state solution}
  \item{Solved }{a logical, TRUE when steady-state has been reached }
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{N      <- 100
Depth  <- seq(0.05, by = 0.1, len = 100)
out    <- OMEXDIAsteady()

# Steady-state concentrations in sediment
CONC  <- out$steady

FDET  <- CONC[1:N]
SDET  <- CONC[(N+1)  :(2*N)]
O2    <- CONC[(2*N+1):(3*N)]
NO3   <- CONC[(3*N+1):(4*N)]
NH3   <- CONC[(4*N+1):(5*N)]
ODU   <- CONC[(5*N+1):(6*N)]

TOC  <- (FDET+SDET)*1200/10^9/2.5     # % organic carbon (excess)

par(mfrow=c(2, 2))
plot(TOC, Depth, ylim = c(10, 0), xlab = "procent", main = "TOC",
        type = "l", lwd=2)
plot(O2, Depth, ylim = c(10, 0), xlab = "mmol/m3", main = "O2",
        type = "l", lwd = 2)
plot(NO3, Depth, ylim = c(10, 0), xlab = "mmol/m3", main = "NO3",
        type = "l", lwd = 2)
plot(NH3, Depth, ylim = c(10, 0), xlab = "mmol/m3", main = "NH3",
        type = "l", lwd = 2)


mtext(outer=TRUE,side=3,line=-2,cex=1.5,"OMEXDIAmodel")}
\references{Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
A model of early diagenetic processes from the shelf to abyssal depths.
Geochimica Cosmochimica Acta, 60(6):1019-1040.

Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
Dynamic response of deep-sea sediments to seasonal variation: a model.
Limnol. Oceanogr. 41(8): 1651-1668.}
\details{The model application just estimates the steady-state condition of the OMEXDIA model, the parameter values are such that there
is almost no overlap between the oxic and anoxic zone of the sediment.

For efficiency reasons, the OMEXDIA diagenetic model was written in Fortran, and this code linked to
the package.}
\keyword{misc}

