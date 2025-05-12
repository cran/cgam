\name{plotpersp_control}
\alias{plotpersp_control}
\title{Control Parameters for plotpersp.cgam}
\description{
Constructs a list of control parameters to customize the appearance of perspective plots created by \code{plotpersp.cgam}.
}
\usage{
plotpersp_control(surface = "mu", x1nm = NULL, x2nm = NULL, categ = NULL,
  col = NULL, random = FALSE, ngrid = 12, xlim = NULL, ylim = NULL, zlim = NULL,
  xlab = NULL, ylab = NULL, zlab = NULL, th = NULL, ltheta = NULL, 
  main = NULL, sub = NULL, ticktype = "simple")
}
\arguments{
  \item{surface}{Character string indicating the surface type to plot. Options are \code{"mu"} for fitted means or \code{"eta"} for systematic component.}
  \item{x1nm}{Optional character strings specifying names of the predictors for the x axis.}
  \item{x2nm}{Optional character strings specifying names of the predictors for the y axis.}
  \item{categ}{Optional character string naming a categorical covariate in the \code{cgam} fit to stratify surfaces. Default is \code{NULL}.}
  \item{col}{Colors used for the surface(s). If \code{NULL}, defaults are applied automatically depending on the number of surfaces.}
  \item{random}{Logical; if \code{TRUE}, colors are randomized for multiple surfaces. Default is \code{FALSE}.}
  \item{ngrid}{Integer; number of grid points along each axis. Default is \code{12}.}
  \item{xlim, ylim, zlim}{Optional numeric vectors of length 2 specifying limits for the x, y, and z axes.}
  \item{xlab, ylab, zlab}{Axis labels for x, y, and z axes. Default is based on the predictor names and surface type.}
  \item{th}{Angle defining the azimuthal direction for viewing the plot (theta).}
  \item{ltheta}{Angle defining the direction of lighting (light theta).}
  \item{main}{Main title of the plot.}
  \item{sub}{Subtitle of the plot.}
  \item{ticktype}{Character string indicating type of tick marks: \code{"simple"} (default) or \code{"detailed"}.}
  %\item{...}{Extra arguments to \code{\link[graphics]{persp}.}
}
\value{
A named list of control settings for use with \code{\link{plotpersp}}.
}
\examples{
ctrl <- plotpersp_control(col = "topo", ngrid = 20, th = 45)
# Then used inside plotpersp:
# plotpersp(fit, control = ctrl)
}
\seealso{
\code{\link{plotpersp}}, \code{\link[graphics]{persp}}
}
\author{
Mary C. Meyer and Xiyue Liao
}
\keyword{graphics}
\keyword{3D}
