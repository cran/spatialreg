\name{aple.plot}
\alias{aple.plot}
\alias{localAple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate profile-likelihood estimator (APLE) scatterplot}
\description{
  A scatterplot decomposition of the approximate profile-likelihood estimator, and a local APLE based on the list of vectors returned by the scatterplot function.
}
\usage{
aple.plot(x, listw, override_similarity_check=FALSE, useTrace=TRUE, do.plot=TRUE, ...)
localAple(x, listw, override_similarity_check=FALSE, useTrace=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a zero-mean detrended continuous variable}
  \item{listw}{a \code{listw} object from for example \code{spdep::nb2listw}}
  \item{override_similarity_check}{default FALSE, if TRUE - typically for row-standardised weights with asymmetric underlying general weights - similarity is not checked}
  \item{useTrace}{default TRUE, use trace of sparse matrix \code{W \%*\% W} (Li et al. (2010)), if FALSE, use crossproduct of eigenvalues of \code{W} as in Li et al. (2007)}
  \item{do.plot}{default TRUE: should a scatterplot be drawn}
  \item{\dots}{other arguments to be passed to \code{plot}}
}
\details{
  The function solves a secondary eigenproblem of size n internally, so constructing the values for the scatterplot is quite compute and memory intensive, and is not suitable for very large n.
}
\value{
\code{aple.plot} returns list with components:
  \item{X}{A vector as described in Li et al. (2007), p. 366.}
  \item{Y}{A vector as described in Li et al. (2007), p. 367.}

\code{localAple} returns a vector of local APLE values.
}
\references{Li, H, Calder, C. A. and Cressie N. A. C. (2007) Beyond Moran's I: testing for spatial dependence based on the spatial autoregressive model. Geographical Analysis 39, pp. 357-375; Li, H, Calder, C. A. and Cressie N. A. C. (2012) One-step estimation of spatial dependence parameters: Properties and extensions of the APLE statistic, Journal of Multivariate Analysis 105, 68-84.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{aple}}}
\examples{
\dontrun{
wheat <- st_read(system.file("shapes/wheat.gpkg", package="spData")[1], quiet=TRUE)
nbr1 <- spdep::poly2nb(wheat, queen=FALSE)
nbrl <- spdep::nblag(nbr1, 2)
nbr12 <- spdep::nblag_cumul(nbrl)
cms0 <- with(as.data.frame(wheat), tapply(yield, c, median))
cms1 <- c(model.matrix(~ factor(c) -1, data=wheat) \%*\% cms0)
wheat$yield_detrend <- wheat$yield - cms1
plt_out <- aple.plot(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
 spdep::nb2listw(nbr12, style="W"), cex=0.6)
lm_obj <- lm(Y ~ X, plt_out)
abline(lm_obj)
abline(v=0, h=0, lty=2)
zz <- summary(influence.measures(lm_obj))
infl <- as.integer(rownames(zz))
points(plt_out$X[infl], plt_out$Y[infl], pch=3, cex=0.6, col="red")
crossprod(plt_out$Y, plt_out$X)/crossprod(plt_out$X)
wheat$localAple <- localAple(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
 spdep::nb2listw(nbr12, style="W"))
mean(wheat$localAple)
hist(wheat$localAple)
opar <- par(no.readonly=TRUE)
plot(wheat[,"localAple"], reset=FALSE)
text(st_coordinates(st_centroid(st_geometry(wheat)))[infl,], labels=rep("*", length(infl)))
par(opar)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

