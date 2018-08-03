#
# Test insol package fDL vs. fDLfun
#

require(insol)
require(GPPFourier)

lat <- 66 # 51.176
lon <- 4.326

insolDL <- NULL
fDLfunDL <- NULL
days <- as.POSIXct("2018-01-01")+(1:365 - 1)*86400
for (day in 1:365){
  # insolDL <- c(insolDL, daylength(lat, lon, unclass(julian(days[day])), 2)[,"daylen"])
  insolDL <- c(insolDL, daylength(lat, lon, JDymd(as.POSIXlt(days[day])$year+1900,as.POSIXlt(days[day])$mon+1, as.POSIXlt(days[day])$mday), 2)[,"daylen"])
  fDLfunDL <- c(fDLfunDL, fDLfun(days[day], phi=lat, lambda=lon))
}

plot(days, insolDL)
lines(days, fDLfunDL*24, col="red")

insolDL[round(days, units="days")==as.POSIXct("2018-08-02")]
fDLfunDL[round(days, units="days")==as.POSIXct("2018-08-02")]*24

insolDL[round(days, units="days")==as.POSIXct("2018-02-01")]
fDLfunDL[round(days, units="days")==as.POSIXct("2018-02-01")]*24