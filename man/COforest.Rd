\name{COforest}
\alias{COforest}
\docType{data}
\title{
Colorado Forest Data Set 
}
\description{
This data set contains 9167 records of different species of live trees for 345 sampled forested plots measured in 2015. 
}
\usage{data("COforest")}
\format{
  A data frame with 9167 observations on the following 19 variables.
  \describe{
    \item{\code{PLT_CN}}{Unique identifier of plot}
    \item{\code{STATECD}}{State code using Bureau of Census Federal Information Processing Standards (FIPS)}
    \item{\code{COUNTYCD}}{County code (FIPS)}
    \item{\code{ELEV_PUBLIC}}{Elevation (ft) extracted spatially using LON_PUBLIC/LAT_PUBLIC}
    \item{\code{LON_PUBLIC}}{Fuzzed longitude in decimal degrees using NAD83 datum}
    \item{\code{LAT_PUBLIC}}{Fuzzed latitude in decimal degrees using NAD83 datum}
    \item{\code{ASPECT}}{a numeric vector}
    \item{\code{SLOPE}}{a numeric vector}
    \item{\code{SUBP}}{Subplot number}
    \item{\code{TREE}}{Tree number within subplot}
    \item{\code{STATUSCD}}{Tree status (0:no status; 1:live tree; 2:dead tree; 3:removed)}
    \item{\code{SPCD}}{Species code}
    \item{\code{DIA}}{Current diameter (in)}
    \item{\code{HT}}{Total height (ft): estimated when broken or missing}
    \item{\code{ACTUALHT}}{Actual height (ft): measured standing or down}
    \item{\code{HTCD}}{Height method code (1:measured; 2:actual measured-length estimated; 3:actual and length estimated; 4:modeled}
    \item{\code{TREECLCD}}{Tree class code (2:growing-stock; 3:rough cull; 4:rotten cull)}
    \item{\code{CR}}{Compacted crown ratio (percentage)}
    \item{\code{CCLCD}}{Crown class (1:open grown; 2:domimant; 3:co-dominant; 4:intermediate; 5:overtopped)}
  }
}
\details{
More detailed explanation can be found in \url{https://www.fia.fs.fed.us/library/database-documentation/current/ver611/FIADB_User_Guide_P2_6-1-1_final.pdf}.
}
\references{
X. Liao and M. Meyer (2019). Estimation and Inference in Mixed-Effect Regression Models using Shape Constraints, with Application to Tree Height Estimation. (\emph{to appear in Journal of the Royal Statistical Society. Series C: Applied Statistics})
}
\source{
It is provided by Forest Inventory Analysis (FIA) National Program.
}
\examples{
\dontrun{
  library(dplyr)
  library(tidyr)
  data(COforest)

  #re-grouping classes of CCLCD: 
  #combine dominant (2) and co-dominant (3) 
  #combine intermediate (4) and overtopped (5)
  COforest = COforest %>% mutate(CCLCD = replace(CCLCD, CCLCD == 2, 3)) %>% 
    mutate(CCLCD = replace(CCLCD, CCLCD == 5, 4))

  #make a list of species, each element is a small data frame for one species
  species = COforest %>% group_by(SPCD) %>% nest()

  #get the subset for quaking aspen, which is the 4th element in the species list
  sub = species$data[[4]]
  #for quaking aspen, there are only two crown classes: dominant/co-dominant 
  #and intermediate/overtopped
  table(sub$CCLCD)
  #   3    4
  #1400  217
  #for quaking aspen, there are only two tree clases: growing-stock and rough cull
  table(sub$TREECLCD)
  #2    3
  #1591   26

  #fit the model
  ansc = cgamm(log(HT)~s.incr.conc(DIA)+factor(CCLCD)+factor(TREECLCD)
  +(1|PLT_CN), reml=TRUE, data=sub)

  #check which classes are significant
  summary(ansc)

  #fixed-effect 95% confidence interval
  newData = data.frame(DIA=sub$DIA,CCLCD=sub$CCLCD,TREECLCD=sub$TREECLCD)
  pfit = predict(ansc, newData,interval='confidence')
  lower = pfit$lower
  upper = pfit$upper
  #we need to use exp(muhat) later in the plot
  muhat = pfit$fit

  #get x and y
  x = sub$DIA
  y = sub$HT

  #get TREECLCD and CCLCD
  z1 = sub$TREECLCD
  z2 = sub$CCLCD

  #plot fixed-effect confidence intervals
  plot(x, y, xlab='DIA (m)', ylab='HT (m)', ylim=c(min(y),max(exp(upper))+10),type='n')
  lines(sort(x[z2==3&z1==2]), (exp(pfit$fit)[z2==3&z1==2])[order(x[z2==3&z1==2])], 
  col='slategrey', lty=1, lwd=2)
  lines(sort(x[z2==3&z1==2]), (exp(pfit$lower)[z2==3&z1==2])[order(x[z2==3&z1==2])], 
  col='slategrey', lty=1, lwd=2)
  lines(sort(x[z2==3&z1==2]), (exp(pfit$upper)[z2==3&z1==2])[order(x[z2==3&z1==2])], 
  col='slategrey', lty=1, lwd=2)

  lines(sort(x[z2==4&z1==2]), (exp(pfit$fit)[z2==4&z1==2])[order(x[z2==4&z1==2])], 
  col="blueviolet", lty=2, lwd=2)
  lines(sort(x[z2==4&z1==2]), (exp(pfit$lower)[z2==4&z1==2])[order(x[z2==4&z1==2])], 
  col="blueviolet", lty=2, lwd=2)
  lines(sort(x[Cz2==4&z1==2]), (exp(pfit$upper)[Cz2==4&z1==2])[order(x[Cz2==4&z1==2])], 
  col="blueviolet", lty=2, lwd=2)

  lines(sort(x[Cz2==3&z1==3]), (exp(pfit$fit)[Cz2==3&z1==3])[order(x[Cz2==3&z1==3])], 
  col=3, lty=3, lwd=2)
  lines(sort(x[Cz2==3&z1==3]), (exp(pfit$lower)[Cz2==3&z1==3])[order(x[Cz2==3&z1==3])], 
  col=3, lty=3, lwd=2)
  lines(sort(x[Cz2==3&z1==3]), (exp(pfit$upper)[Cz2==3&z1==3])[order(x[Cz2==3&z1==3])], 
  col=3, lty=3, lwd=2)

  lines(sort(x[Cz2==4&z1==3]), (exp(pfit$fit)[Cz2==4&z1==3])[order(x[Cz2==4&z1==3])], 
  col=2, lty=4, lwd=2)
  lines(sort(x[Cz2==4&z1==3]), (exp(pfit$lower)[Cz2==4&z1==3])[order(x[Cz2==4&z1==3])], 
  col=2, lty=4, lwd=2)
  lines(sort(x[Cz2==4&z1==3]), (exp(pfit$upper)[Cz2==4&z1==3])[order(x[Cz2==4&z1==3])], 
  col=2, lty=4, lwd=2)
  legend('bottomright', bty='n', cex=.9,c('dominant/co-dominant and growing stock',
  'intermediate/overtopped and growing stock','dominant/co-dominant and rough cull',
  'intermediate/overtopped and rough cull'), col=c('slategrey',"blueviolet",3,2), 
  lty=c(1,2,3,4),lwd=c(2,2,2,2), pch=c(24,23,22,21))
  title('Quaking Aspen fits with 95% CI', cex=1.5)
}


}
\keyword{datasets}
