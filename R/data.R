#' White Crown Sparrow Song Observations
#' 
#' A dataset of mean trill bandwidth from 3 WCSP male populations in 1969 and 2005 
#' 
#' @format  A data frame with 89 rows and 3 variables:
#' \describe{
#'   \item{Population}{Locality of the observed male sparrow}
#'   \item{Year}{Year in which samples were taken}
#'   \item{Trill.FBW}{Mean observed trill bandwidth}
#' }
#' @source \url{http://FIXME/}
"song.data"

#' Default model parameters
#' 
#' Default model parameters for WCSP male trill bandwidth in 1969
#' 
#' @format  A data frame with 89 rows and 3 variables:
#' \describe{
#'   \item{learning.error.d}{Drift in trill bw due to learning error}
#'   \item{learning.error.sd}{Dispersal in trill bw due to learning error}
#'   \item{global.parms$learning.error.d}{NULL}
#'   \item{n.territories}{Number of available territories}
#'   \item{mortality.a}{Mortality a}
#'   \item{mortality.j}{Mortality J}
#'   \item{lifespan}{Mean WCSP lifespan}
#'   \item{phys.lim.min}{Minimum possible trill bandwidth}
#'   \item{phys.lim.max}{Maximum possible trill bandwidth}
#'   \item{male.fledge.n.mean}{Mean number of males fledged}
#'   \item{male.fledge.n.sd}{Std deviation of the number of males fledged}
#'   \item{disp.age}{Age at which males disperse}
#'   \item{disp.distance.mean}{Mean dispersal distance}
#'   \item{disp.distance.sd}{Std deviation of dispersal distance}
#'   \item{terr.turnover}{Relative number males who lose their territory}
#'   \item{male.fledge.n}{vector of length n.territories with number of males fledged for each territory}
#' }
"glo.parms"