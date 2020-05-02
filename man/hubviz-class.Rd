
\name{hubviz-class}
\alias{hubviz-class}
\docType{class}

\title{Class "hubviz"}
\description{
This class represents hubviz fit.
}


\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hubviz", ...)}.
}

\section{Slots}{
\describe{

\item{\code{data}}{Object of class \code{"matrix"}, 
    	representing the raw data for hubViz model. }

\item{\code{init}}{Object of class \code{"list"}, 
    	representing the setting for hubViz model fitting.} 

\item{\code{result}}{Object of class \code{"list"}, 
    	representing the fitted hubViz model.} 
}}


\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "hubviz")}: provide brief summary of the object. }
    
    \item{estimate}{\code{signature(object = "hubviz")}: 
		provide estimate of latent posisions.}
		
    \item{plot}{\code{signature(x = "hubviz", y=missing)}: 
    provide hub-centric visualization}

    \item{hpdplot}{\code{signature(object = "hubviz")}: 
		provide HPD plot for latent posisions.}

}}

\details{

}

\references{
Nam JH, Yun J, Jin IH and Chung D (2020) ''hubViz: A novel tool for hub-centric visualization''.
}

\author{ Jin Hyun Nam, Ick Hoon Jin }

\seealso{
\code{\link{hubviz}}.
}

\examples{
showClass("hubviz")

n <- 100
p <- 10
bmat <- matrix( 0, n, p )
for ( j in 1:ncol(bmat) ) {
 bmat[ 1:(10*j), j ] <- 1
}
rownames( bmat ) <- paste( "Sample", 1:n, sep="" )
colnames( bmat ) <- paste( "V", 1:p, sep="" )

nsamp <- nrow( bmat ); nitem <- ncol( bmat )
fit <- hubviz( bmat, nsample = nsamp, nitem = nitem)
fit

estimate( fit )

plot( fit )

hpdplot( fit )
}

\keyword{classes}