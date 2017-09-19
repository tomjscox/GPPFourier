#' Relative amount of daylight hours at a specified date and location.
#' @param date POSIXct object or date specificied in unambiguous format. See \code{\link{as.POSIXct}}
#' @param phi Latitude
#' @param lambda Longitude
#' @param H Height of location where fDL is to be calculated
#'
#' @return The fraction of daylight hours at the specified date and location. Sunrise and Sunset are calculated with \code{\link{SunRiseSet}}.
#' @export fDLfun
#' @examples fDLfun("2016-06-21")
#' @seealso \code{\link{SunRiseSet}}
#' @author Tom Cox <tom.cox@uantwerp.be>

fDLfun <- function(date="2016-07-01", phi=51.176, lambda=4.326, H=0) {
date <- round.POSIXt(as.POSIXct(date), units="days")
if (length(date)<2){
  DL <- diff(SunRiseSet(date,phi=phi,lambda=lambda))
} else {
  AllDL <- NULL
  dates <- unique(date)
  for (i in 1:length(dates)) {
    DL <- diff(SunRiseSet(dates[i],phi=phi,lambda=lambda))
    AllDL <- c(AllDL, DL)
    }
  DL <- AllDL[match(date, dates)]
  }
fDL <- DL/24
return(fDL)
}



#' Calculate sunset and sunrise for a given date at given location and observer height
#' 
#' @param date POSIXct object or date specificied in unambiguous format. See \code{\link{as.POSIXct}}
#' @param phi Latitude
#' @param lambda Longitude
#' @param H Height of location where fDL is to be calculated
#'
#' @return Two element vector with sunrise and sunset time
#' @details This is an implementation of the algorithm found in the explanatory supplement of the astronomical almanac p484 and further.
#' @examples fDLfun("2016-06-21")
#' @seealso \code{\link{fDLfun}}
#' @author Tom Cox <tom.cox@uantwerp.be>
#' 
#' @export SunRiseSet

SunRiseSet <- function(date="2016-07-01", phi=51.176, lambda = 4.326, H=0){

# Algorithm can be checked with e.g.
# http://www.timeanddate.com/worldclock/astronomy.html?n=48&month=1&year=1997&obj=sun&afl=-12&day=1

# Calculation of Sunrise/Sunset hour
# Formulas from the explanatory supplement of the astronomical almanac p484 and further
# phi is latitude, default values for Kruibeke
# lambda is longitude
# H is height of the observer

# For the sun
# JD = Julian date
# UT = Universal Time

date <- as.POSIXct(date)
date <- trunc(date,units="days")

JulianDay <- unclass(2440587.50000 + (julian((date))))[1] #Julian date of 1-1-1970 taken from http://www.usno.navy.mil/USNO/astronomical-applications/data-services/jul-date
JD <- c(JulianDay, JulianDay)	
#initial guess of UT <- 12h

UT0 <- c(0,0)
UT <- c(12,12)

#iterative solution
#replace UT by a value between 0 and 24 by adding multiples of 24
#replace UT0 by UT until the difference between the two is less than 0h.008

while (Mod(UT[1]-UT0[1])>.008 | (Mod(UT[2]-UT0[2])>.008))
{
UT0 <- UT
T <- (JD + UT0/24 - 2451545.0)/36525	#T is number of centuries from J2000

#Calculate Solar arguments:
L <- (280.460 + 36000.770*T)%%360	#remove multiples of 360 deg. first term and multiplication factor of T are in deg
G <- 357.528 + 35999.050*T
lambda1 <- L + 1.915*sin(2*pi/360*G) + 0.020*sin(2*pi/360*2*G)
eps <- 23.4393 - 0.01300*T
E <- -1.915*sin(2*pi/360*G)-0.020*sin(2*pi/360*2*G)+2.466*sin(2*pi/360*2*lambda1) - 0.053*sin(2*pi/360*4*lambda1)
GHA <- 15*UT - 180 + E
delta <- asin(sin(2*pi/360*eps) * sin(2*pi/360*lambda1))
SD <- 0.267/(1-0.017*cos(2*pi/360*G))


h <- -50/60 - 0.0353*sqrt(H)	#H is height of the observer

cost <- (sin(2*pi/360*h) - sin(2*pi/360*phi)*sin(delta))/(cos(2*pi/360*phi)*cos(delta))
t <- c(0,0)
for (i in 1:2){
if (cost[i] >1) {t[i] <- 0
} else if (cost[i] < -1) {t[i] <- pi
} else t[i] = 360/2/pi*acos(cost[i])
}

UT = UT0 - (GHA + lambda + c(t[1],-t[2]))/15 # +t is for sunrise and -t is for set 

UT <- UT%%24

}
return(UT)
}
