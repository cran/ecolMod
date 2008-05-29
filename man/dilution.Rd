\name{dilution}
\alias{dilution}
\title{Draws dilution culture setup}
\description{Draws the framework of a dilution culture. Used as a template to plot the 
	flow diagrams for dilution-type models, where there is a continuous inflow of medium from 
  a vessel into a well-stirred tank (the culture vessel). The volume in the culture tank stays constant.}
\usage{dilution(main=c("",""),int="")}
\arguments{
  \item{main }{main text, consiting of two strings, positioned above medium vessel (1st string) and culture vessel (2nd string) of dilution diagram plot}
  \item{int }{text above the dripping outlets}
}
\value{
  a list containing:
  \item{p1 }{two-valued vector with the x-y positions of the middle of the (large) culture vessel}
  \item{p2 }{two-valued vector with the x-y positions of the middle of the (small) medium vessel}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{par(mar=c(0,0,0,0))                         
dd<-dilution(main=c("Stock","Stirred tank"),int="Flow,Q")
text(dd$p2[1],dd$p2[2],"Ain",font=2)
text(dd$p1[1],dd$p1[2]+0.03,"Volume V",font=2)
text(dd$p1[1],dd$p1[2]-0.03,"[A]",font=2)}
\keyword{aplot}

