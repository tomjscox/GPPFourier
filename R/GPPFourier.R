#
# Core functions to calculate GPP from O2 time series making use of Fourier method
#

#'
#' @title Pre-process O2 data for GPP calculation
#' @description Removes longer term trend of time series and confine series to integer
#' number of days
#' @param x A vector containing consecutive O2 concentrations, sampled at regular intervals (time step = dt).
#' @param dt Sampling time step. Class {difftime}, or \code{numeric} when units are provided.
#' @param units Units of sampling time step
#' @param Detrend Toggle time series detrending
#' @param filter Toggle time series filtering
#' @param Nfilt Moving average filter width (see documentation of
#' \code{filter()})
#' @param circular Moving average boundary condition (see documentation of
#' \code{filter()})
#' @param sides Moving average central or one sided (see documentation of
#' \code{filter()})
#' @details GPPFourierPreprocess() is called by GPPFourier() to pre-process the time series before calculating the Fourier amplitude at diel frequency and estimate GPP. 
#' @return A list with components \item{filt}{The long term trend}
#' \item{res}{The original series without the long term trend, confined to an
#' integer number of days} \item{indices}{The indices of the subset of the
#' series that are retained, after filtering and confining to an integer number
#' of days.}
#' @seealso \code{\link{GPPFourier}}, \code{\link{GPPFourier_t}}, \code{\link{WindowGPPFourier}},
#' \code{\link{WindowGPPFourier.gts}}
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @references Cox T.J.S. et al. (2015) 'Estimating primary production from oxygen time series: a novel approach in the frequency domain', Limnology And Oceanography:Methods 13, 529-552. DOI: 10.1002/lom3.10046

#' @keywords utilities
#' @export GPPFourierPreprocess
GPPFourierPreprocess <- function(x,		# O2 time series from which GPP is to be calculated
			dt=NULL,		# Sampling time step
			units=c("days", "hours", "mins"),		# Time units of sampling time step
			Detrend=TRUE,		# Toggle time series detrending
			filter = TRUE, 		# Toggle time series filtering
			Nfilt=NULL,		# Moving average filter width (see documentation of filter())
			circular = FALSE, 	# Moving average boundary condition (see documentation of filter())
			sides = 2)		# Moving average central or one sided (see documentation of filter()) 


{

# Removes longer term trends and confines time series to integer number of periods 
# 
#

unit <- match.arg(units)

if (is.null(dt)) stop("Provide dt")
if (!inherits(dt,"difftime")){
if (is.null(units)) {
  stop("Either provide dt as difftime object, or provide units of dt")
} else {
  dt <- as.difftime(dt, units=unit)
}}

if (filter&is.null(Nfilt)) stop("Provide filter length (Nfilt)")
if (Nfilt>length(x)) stop("Filter longer than time series")


dt <- as.double(dt, units="days")

# confine series to integer number of periods (T=1 day = 1/dt samples)

Nperiod <- floor(length(x)*dt)/dt
indices <- 1:Nperiod
x <- x[indices] 

if (Detrend){
detrendx <- detrend(x, returntrend=TRUE)
x <- detrendx$detrended
trendx <- detrendx$trend
} else {
trendx <- rep(mean(x), length(x))
}

if (filter){

cfilt <- rep(1/Nfilt,Nfilt)

filtx <- filter(x,cfilt,sides=sides,circular=circular)
resx <- x - filtx
indices <- indices[!is.na(filtx)]

} else {
  filtx <- rep(mean(x),length(x))
  resx <- x - mean(x)
}

# confine filtered series to integer number of periods (T=1 day = 1/dt samples)
filtx <- filtx[!is.na(filtx)] 
resx <- resx[!is.na(resx)]

Nperiod <- floor(length(filtx)*dt)/dt
indices <- indices[1:Nperiod]
filtx <- filtx[1:Nperiod] 
resx <- resx[1:Nperiod]

return(list(filt = filtx, res=resx, trend=trendx, indices=indices))
}


#'
#' @title Calculate GPP from O2 time series
#' 
#' @description This is the workhorse to calculate aquatic Gross Primary Production from
#' high frequency dissolved oxygen data. It will often be more convenient to use the higher
#' level functions \code{\link{WindowGPPFourier}} and
#' \code{\link{WindowGPPFourier.gts}}
#' 
#' @param x Regularly sampled O2 time series. x can be a dataframe with time (POSIXt) in first column and O2 concentrations in 
#' second column or a vector of concentrations sampled with time step = dt
#' @param dt 	Sampling time step. Not used when x is a data frame, then dt is calculated #' from the time spacing between 
#' the first two samples. Either a difftime object or a numerical value. When dt is given as a numerical value, \code{units} 
#' have to be provided.
#' @param units Unit of sampling time step
#' @param Nfilt Moving average filter width
#' @param fDL Relative fraction of light hours during the day. When \code{x} is a data frame, and no value for fDL is provided, 
#' fDL is calculated with \code{\link{fDLfun}}
#' @param filtcorrect Logical controlling whether GPP estimate
#' is corrected for signal falsely removed by filtering
#' @param ... Other parameters to be passed to \code{\link{GPPFourierPreprocess}}
#' @param phi Optional: Lattitude. Only used when fDL=NULL and daylength is calculated with \code{\link{fDLfun}}
#' @param lambda Optional: Lattitude. Only used when fDL=NULL and daylength is calculated with \code{\link{fDLfun}}
#' @return Volume specific average gross primary production in units of oxgen as in \code{x}. 
#' per time unit (unit of \code{dt}, given by \code{units}). 
#' Specifically: if O2 concentrations are in mg/L, and \code{units = days}, calculated GPP has units mg/L/day.
#' If surface specific primary production (e.g. per m2) is desired, this number has to be multiplied by the depth of the water body.
#' @details no details
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @seealso \code{\link{GPPFourierPreprocess}}, \code{\link{GPPFourier_t}}, \code{\link{WindowGPPFourier}},
#' \code{\link{WindowGPPFourier.gts}}
#' @references Cox T.J.S. et al. (2015) 'Estimating primary production from oxygen time series: a novel approach in the frequency domain', Limnology And Oceanography:Methods 13, 529-552. DOI: 10.1002/lom3.10046
#'
#' @keywords utilities
#' @examples
#' DO <- Kruibeke[Kruibeke$time>="2010-06-03"&Kruibeke$time<="2010-06-13",]
#' dt <- as.numeric(difftime(DO$time[2],DO$time[1],units="days"))
#' DL <- fDLfun(DO$time[1], phi=51.176, lambda=4.326) 
#' GPP <- GPPFourier(DO$O2,dt=dt,Nfilt=100,fDL=DL)
#' 
#' @export GPPFourier

GPPFourier <- function(	x,			# O2 time series from which GPP is to be calculated
			dt=NULL,		# Sampling time step
			units=c("days","hours","mins"),		# Time units of sampling time step
			Nfilt=NULL,		# Moving average filter width
			fDL=NULL,		# Relative fraction of light hours during the day
			filtcorrect = FALSE,	# Logical controlling whether GP estimate is corrected for signal removed by moving average filtering 
			phi=NULL,     # Optional: Lattitude. Only used when fDL=NULL and daylength is calculated with \code{\link{fDLfun}}
			lambda=NULL,  # Optional: Longitude. Only used when fDL=NULL and daylength is calculated with \code{\link{fDLfun}}
			...)			# Other parameters to be passed to GPPFourierprocess
{
# Calculate average gross production from O2 time series
# Units: unit of oxgen as in ts per days
#

#
# Error handling
# 

ts <- x # internal alias

if (inherits(ts,"data.frame")){
  if (any(is.na(ts[,c(1,2)]))) stop("NAs not allowed")
  if (dim(ts)[2]>=2){
    t <- ts[,1]
    if (inherits(t, "POSIXt")){ 
      if (any(0!=diff(diff(round.POSIXt(t, units="secs"))))) stop("Samples are not equidistant and/or gaps are present in the series. Re-sample, fill gaps and/or split series in continuous subsets")    
      dt <- difftime(t[2],t[1],units="days")        
    } else {
      if (any(0!=round(diff(diff(t)), digits=5))) stop("Samples are not equidistant and/or gaps are present in the series. Re-sample, fill gaps and/or split series in continuous subsets")      
      dt <- diff(c(t[1],t[2]))  
    }
    if (is.null(fDL)){
      if (is.null(phi)|is.null(lambda)) stop("Provide longitude and lattitude") else {
      fDL <- fDLfun(t, phi=phi, lambda=lambda)
      fDL <- mean(fDL)    # average fDL of time series is used
      }
    }
    
    O2 <- ts[,2]
  }
} else {
  if (is.null(dt)) stop("Provide dt")
  if (is.null(fDL)) stop("Provide daylight fraction of day (fDL)")
  if (any(is.na(ts))) stop("NAs not allowed")
  O2 <- ts
  # fDL <- rep(fDL, length(O2))
  unit <- match.arg(units)
  if (!inherits(dt,"difftime")){
    if (is.null(units)) {
      stop("Either provide dt as difftime object, or provide units of dt")
    } else {
      dt <- as.difftime(dt, units=unit)
    }
  }
  t <- 1:length(O2)*dt
}
  

# if (!inherits(dt,"difftime")){
# if (is.null(units)) {
#   stop("Either provide dt as difftime object, or provide units of dt")
# } else {
#   dt <- as.difftime(dt, units=unit)
# }}

# if (is.null(Nfilt)) stop("Provide moving average filter width (Nfilt)")
# if (is.null(fDL)) stop("Provide daylight fraction of day (fDL)")
# if (as.numeric(length(x) < Nfilt)) stop("Filter longer than time series")


# angular frequency (rad per day) associated with diurnal periodicity, i.e. T = 1 day
w <- 2*pi/1

# Remove longer term trends and confines time series to integer number of periods 
xpreprocessed <- GPPFourierPreprocess(x, dt=dt, Nfilt=Nfilt, ...)

filtx <- xpreprocessed$filt
resx <- xpreprocessed$res
indices <- xpreprocessed$indices

dt <- as.double(dt,units="days")		
# direct calculation of A24
times <- (1:length(resx))*dt
A24 <- 2*sum(exp(-times*1i*w)*resx)/length(resx)

# correction for truncated sinusoidal form of GPP
theta1 <- pi*(1-fDL)
truncFac <- (sin(theta1) + (pi-theta1) * cos(theta1))/(pi-theta1 + 0.5 * sin(2*theta1))
A24 <- A24*truncFac

# Calculate time derivative in frequency domain
A24 <- 1i*w*A24

# Correct for frequency response of moving average filter
if (filtcorrect){
if (is.null(Nfilt)){
warning("Nfilt=NULL. No filter correction calculated")
} else { 
L <- Nfilt
maxfreq <- (length(resx) - 1)/(length(resx))/dt
w24h <- 2*pi/maxfreq
A24 <- A24/(1-Hw(w24h,Nfilt))
}}

return(Mod(A24))
}



#' 
#' @title Calculate GPP from O2 time series with gaps, in consecutive time-windows of prescribed
#' length
#' 
#' @description Calculate GPP from O2 time series in consecutive time-windows of prescribed
#' length. High level function to handle time series with gaps. Short gaps are
#' interpolated, series is split at large gaps. 
#' 
#' 
#' @param x Dataframe containing time (POSIXt) in first column and O2 concentrations in the second column. 
#' @param dt Sampling time step. Either a difftime object or a numerical value. 
#' When \code{dt} is given as a numerical value, \code{units} have to be provided (defaults to \code{days}). When omitted \code{dt} is calculated from the time spacing between the first two samples in the data frame.
#' @param units Unit of sampling time step
#' @param gapMaxN Minimum number of missing data points to be considered a gap in the time series. Gaps smaller than gapMaxN are interpolated. Defaults to the number of samples corresponding to a 4 hour gap.
#' @param Width Width [days] of the time-windows for which GPP is calculated
#' @param filtWidth Length of moving average filter [hours] to filter O2
#' @param ... Other parameters to be passed to WindowGPPFourier()
#' @return Average gross primary production in units of oxgen as in \code{x}
#' per day.
#' @importFrom Tides gapsts
#' @details This high level function analyses multiple regular O2 time series at once, each separated by a gap larger than gapMaxN (the number of . This is convenient e.g. for regular O2 data from a single location, where some data is missing due to sensor maintenance or replacement. \code{x} is assumed to contain all these data, pasted together in a single data frame (e.g. by \code{\link{cbind}}ing all individual regular time series or just reading from a single data file). Additionally, missing values within the regular series are allowed and are interpolated with \code{\link{gapfill}}. Uses \code{\link{WindowGPPFourier}}.
#' @seealso \code{\link{GPPFourierPreprocess}}, \code{\link{GPPFourier_t}}, \code{\link{SunRiseSet}},
#' \code{\link{WindowGPPFourier}}
#' @keywords utilities
#' @references Cox T.J.S. et al. (2015) 'Estimating primary production from oxygen time series: a novel approach in the frequency domain', Limnology And Oceanography:Methods 13, 529-552. DOI: 10.1002/lom3.10046
#' @examples 
#' par(mfrow=c(2,1))
#' plot(Kruibeke, pch=20,xlim=as.POSIXct(c("2010-01-01", "2010-12-31")))
#' GPPAll_4 <- WindowGPPFourier.gts(Kruibeke, 
#'                                  gapMaxN = 10, 
#'                                  Width = 10, 
#'                                  filtWidth=1*24, 
#'                                  phi=51.176,
#'                                  lambda=4.326, 
#'                                  Detrend=TRUE, 
#'                                  filter=TRUE, 
#'                                  filtcorrect=FALSE)
#' plot(GPPAll_4$time,GPPAll_4$GPP*9.3, 
#'      col="black", 
#'     pch=19, 
#'     type="b", 
#'     xlab="time", 
#'     ylab="GPP", 
#'     xlim=as.POSIXct(c("2010-01-01", "2010-12-31")))


#' 
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @export WindowGPPFourier.gts
WindowGPPFourier.gts <- function(x, 
                                 dt = difftime(x[2,1],x[1,1],units="days"),
                                 units = c("days", "hours", "mins"),
                                 gapMaxN = 4/24/dt, 
                                 Width = 14, 
                                 filtWidth=16,
                                 
                                 ...) {
  
  # x: data frame with first column t, and second column O2
  # gapMaxN: number of missing samples that will be interpolated. More than gapMaxN missing samples will be handled as a gap
  # with dt = 10 min, N = 6 amounts to 1h-gaps being linearly interpolated, 6*4 = 4h gaps interpolated
  
  ts <- x		# internal alias
  ts$O2interpol <- ts[,2]
  
  unit <- match.arg(units)
  if (!inherits(dt,"difftime")){
    if (is.null(units)) {
      stop("Either provide dt as difftime object, or provide units of dt")
    } else {
      dt <- as.difftime(dt, units=unit)
    }
  }
  
  dt <- as.double(dt, units="days")
  
  
  
  # Remove all gaps > dt*gapMaxN
  gaps2 <- gapsts(ts[,1], dtMax=dt*gapMaxN, unit="days")
  ts$N <- 0
  ts$N[is.element(ts[,1],gaps2$t2)] <- 1
  ts$N <- cumsum(ts$N)
  
  # interpolate all gaps < dt*gapsMaxN
  tslist <- split(ts[,1:2], ts$N)
  tsfilledAll <- lapply(tslist, gapfill)
  tsfilled <- do.call(rbind, tsfilledAll)
  
  # re-add N
  gaps2 <- gapsts(tsfilled[,1], dtMax=dt*gapMaxN, unit="days")
  tsfilled$N <- 0
  tsfilled$N[is.element(tsfilled$t,gaps2$t2)] <- 1
  tsfilled$N <- cumsum(tsfilled$N)
  
  # Windowed GP calculation on 'continuous' subsets of 10d or more, with small gaps filled
  tslist <- split(tsfilled[,c("t","y")], tsfilled$N)
  dtlist <- lapply(tslist, function(x) {x[1,1]- x[2,1]})
  sublengths <- tapply(tsfilled$t, tsfilled$N, length)
  
  tslist <- tslist[sublengths > Width/dt]
  GPAll <- lapply(tslist, WindowGPPFourier, Width=Width, filtWidth=filtWidth, ...)
  
  GPAll <- do.call(rbind, GPAll)
  return(GPAll)
}



#' 
#' @title Calculate GPP from O2 time series in consecutive time-windows of prescribed
#' length
#' 
#' @description Calculate GPP from O2 time series in consecutive time-windows of prescribed
#' length. 
#' 
#' @param x Vector or dataframe containing consecutive O2 concentrations, sampled at regular intervals (time step = dt). In case of irregular time intervals, consider \code{\link{resample}}.
#' When x is a data frame, time (POSIXt) must be in first column and O2 in second column; when x is a vector, dt must be provided. 
#' @param dt Sampling time step. Either a difftime object or a numerical value. When \code{dt} is given as a numerical value, 
#' \code{units} have to be provided. Can be omitted when \code{x} is a data frame, then \code{dt} is calculated from the 
#' time spacing between the first two samples. 
#' @param units Unit of sampling time step
#' @param Width Width [days] of the time-windows for which GPP is calculated
#' @param Nblocks Number of consecutive subsections on which to calculate GPP. By default the number of windows of Width fitting into the range of x[,1]
#' Default the largest number of blocks ow width=Width that
#' @param phi Latitude of location where O2 series was recorded. Used to
#' calculate relative fraction of light hours during the day, during which
#' production takes place
#' @param lambda Longitude of location where O2 series was recorded
#' @param filtWidth [hours] Lenght of moving average filter to filter O2
#' series. See GPPFourierPreprocess()
#' @param fDL Optional. Relative fraction of daylight. When \code{x} is a data frame, fDL is calculated with \code{\link{fDLfun}} 
#' @param ... Other parameters to be passed to \code{GPPFourier()}
#' @return Average gross primary production in units of oxgen as in \code{x}
#' per day.
#' @details no details
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @seealso \code{\link{GPPFourierPreprocess}}, \code{\link{GPPFourier}}, \code{\link{GPPFourier_t}}, \code{\link{WindowGPPFourier.gts}}, \code{\link{SunRiseSet}}
#' @keywords utilities
#' 
#' @references Cox T.J.S. et al. (2015) 'Estimating primary production from oxygen time series: a novel approach in the frequency domain', Limnology And Oceanography:Methods 13, 529-552. DOI: 10.1002/lom3.10046
#' 
#' 
#' @export WindowGPPFourier
WindowGPPFourier <- function(x, 	 
                             dt=NULL,
                             units=c("days","hours","mins"),	
                             Width=14, 	
                             Nblocks=floor(unclass(difftime(range(ts[,1])[2], range(ts[,1])[1],units="days"))[1]/Width), 
                             fDL=NULL,
                             phi, 		
                             lambda, 
                             filtWidth=16, 	
                             ...){		
  
  ts <- x	# internal alias
  
  if (inherits(ts,"data.frame")){
    if (any(is.na(ts[,c(1,2)]))) stop("NAs not allowed")
    if (dim(ts)[2]>=2){
      t <- ts[,1]
      if (any(0!=diff(diff(round.POSIXt(t, units="secs"))))) stop("Samples are not equidistant and/or gaps are present in the series. Re-sample, fill gaps and/or split series in continuous subsets")    
      dt <- difftime(t[2],t[1],units="days")
      O2 <- ts[,2]
      if (is.null(fDL)){
        fDL <- fDLfun(t, phi=phi, lambda=lambda)
      } else if (length(fDL)==1){
        fDL <- rep(fDL, length(ts[,2]))
      } else if (length(fDL)!= length(ts[,2])){
        stop("When fDL is provided, it should either be single valued or a vector of same length as the time series")
      }
    }
  } else {
    if (is.null(dt)) stop("Provide dt")
    if (any(is.na(ts))) stop("NAs not allowed")
    if (is.null(fDL)) {
      stop("Provide daylight fraction of day (fDL)")
    }    else if (length(fDL)==1){
      fDL <- rep(fDL, length(ts[,2]))
    } else if (length(fDL)!= length(ts[,2])){
      stop("When fDL is provided, it should either be single valued or a vector of same length as the time series")
    }
    
    O2 <- ts
    unit <- match.arg(units)
    if (!inherits(dt,"difftime")){
      if (is.null(units)) {
        stop("Either provide dt as difftime object, or provide units of dt")
      } else {
        dt <- as.difftime(dt, units=unit)
      }
    }
    t <- 1:length(O2)*dt
  }
  
  dtnumeric <- as.double(dt, units="days")
  Nfilt <- filtWidth/dtnumeric/24
  
  
  blockst1 <- t[1] + (0:(Nblocks-1))*86400*Width
  blockmids <- t[1] + Width/2*86400 + (0:(Nblocks-1))*86400*Width
  blocknum <- approx(x=blockst1,y=0:(Nblocks-1),xout=t,method="constant", rule=2)$y
  blocklength <- tapply(t, blocknum, length)
  
  # calculate GPP on all consecutive subsets of Width
  GPAll <- NULL
  for (i in unique(blocknum))
  {
    O2i <- O2[blocknum==i]
    fDLi <- mean(fDL[blocknum==i])
    ti <- t[blocknum==i]
    if (length(O2i) < 1.5*Nfilt){
      GPAll <- c(GPAll,NA)
    } else {
      # RiseSet <- SunRiseSet(ti,phi=phi,lambda=lambda) + 2   #summertime = GMT+2
      # DL <- diff(RiseSet)
      # fDL <- DL/24
      GPO2 <- GPPFourier(x=O2i,dt=dt,Nfilt= Nfilt, fDL=fDLi, ...) 		# Unit: mmol O2/m3/d
      GPAll <- c(GPAll,GPO2)
    }
  }
  return(data.frame(time=blockmids[unique(blocknum+1)],GPP=GPAll, t1=blockst1[unique(blocknum+1)], t2=blockst1[unique(blocknum+1)] + Width*86400))
}



#'
#' @title Complex demodulation
#' @description Retreive the time varying amplitude of a signal with fixed carrier frequency by complex demodulation
#' @param x Vector containing time series observations at constant sampling interval
#' @param dt Time step of time series. Arbitrary units. Class \code{Numeric} or \code{difftime}
#' @param w1 Carrier frequency. Inverse unit of \code{dt}.
#' @param filttype Filter type to be applied for low pass filtering. Only "MA" (moving average) is implemented.
#' @param Nf Moving average filter width
#' @param nf Number of iterative applications of moving average filter. To remove filter bandwidth.
#' @param freq Cut-off frequency of Buttersworth filter
#' @details no details
#' @export demod

demod <- function(	x, 				# Time series to be demodulated
                   dt,				# Time series time step
                   w1 = 2*pi/1, 			# Carrier frequency
                   filttype=c("MA"),		# Filter type
                   Nf, 				# Moving average filter width
                   nf = 1				# Number of moving averaging filter passes
                   )
  
  
{
  
  filttype <- match.arg(filttype)
  
  if (Nf>length(x)) {
    warning("Demodulation filter longer than time series")
    return(rep(NA, length(x)))
  }

  if (inherits(dt, "difftime")) dt <- unclass(dt)
  
  times <- (1:length(x))*dt
  demodfac <- exp(-1i*w1*times)
  xdemod <- (x - mean(x))*demodfac
  
  if (filttype=="MA"){
    demodfiltRe <- Re(xdemod)
    demodfiltIm <- Im(xdemod)
    for (i in 1:nf){
      demodfiltRe <- filter(demodfiltRe,rep(1/Nf,Nf))
      demodfiltIm <- filter(demodfiltIm,rep(1/Nf,Nf))
    }
    demodfilt <- demodfiltRe + demodfiltIm*1i
  } else {
  }
  return(2*demodfilt)
}

#'
#' @title Calculate GPP(t) from by complex demodulation of O2 time series 
#' @description Calculate GPP(t) from by complex demodulation of O2 time series 
#' @param x Regularly sampled O2 time series. x can be a dataframe with time (POSIXt) in first column and O2 concentrations in second column
#' or a vector of concentrations sampled with time step = dt
#' @param dt Sampling time step. Either a difftime object or a numerical value. 
#' When \code{dt} is given as a numerical value, the unit is assumed to be
#' days, unless \code{units} is provided. Can be omitted when \code{x} is a data frame, then \code{dt} is calculated from the time spacing between the first two samples. 
#' @param units Unit of sampling time step. 
#' @param Nfilt Moving average filter width for detrending time series. See \code{\link{GPPFourierPreprocess}}
#' @param NLowPass Moving average filter width for low pass filtering demodulated O2 series. See \code{\link{demod}}
#' @param fDL Optional. Relative fraction of light hours during the day. When \code{x} is a data frame, and no value for fDL is 
#' provided, fDL is calculated with \code{\link{fDLfun}} 
#' @param MAcorrect correction factor for moving average filter
#' @param nf Number of iterative applications of moving average filter, to reduce filter width
#' @param phi Latitude
#' @param lambda Longitude
#' @param ... Other parameters to be passed to GPPFourierPreprocess
#' @seealso \code{\link{GPPFourierPreprocess}}, \code{\link{GPPFourier}}, \code{\link{WindowGPPFourier.gts}}, \code{\link{SunRiseSet}}

#' @return Data frame with time and GPPt column.
# #' @details no details
#' @examples 

#' # Calculate GPP(t) by complex demodulation of simulated water column time series
#' # More examples in vignette("GPPFourier").
#' 
#' dt <- as.numeric(difftime(watercolumn$time[2], watercolumn$time[1] , units="days"))
#' Nfilt <- 1/dt
#' GPPt <- GPPFourier_t(watercolumn[,c("time", "O2")], 
#'                         dt=dt, 
#'                         units="days", 
#'                         Detrend=TRUE, 
#'                         filter=TRUE, 
#'                         Nfilt=Nfilt, 
#'                         NLowPass=Nfilt, 
#'                         fDL=NULL, 
#'                         circular=FALSE, 
#'                         sides=2, 
#'                         nf=1)
#' 
#' 
#' par(mfrow=c(2,1), cex=1.2)
#' plot(watercolumn[,c("time", "O2")], type="l", xlab="", ylab=expression(paste(O[2], " [", mu, "M]")))
#' title(main="Water column")
#' plot(GPPt, type="l", lwd=3, ylim=c(0,30), xlab="", ylab="GPP")
#' lines(watercolumn[,c("time","GPP")], col="red")
#' legend("topleft", 
#'          lty=1, 
#'          col=c("red", "black"), 
#'          legend=c( "Simulated GPP", expression(paste("Complex demodulated ", O[2], " series"))), 
#'          bty="n")
#' 
#' 


#' @references Cox et al (2017). Tune in on 11.57 muHz and listen to primary production. Biogeosciences Discussions. doi:10.5194/bg-2017-81
#' @export GPPFourier_t

GPPFourier_t <- function(	x,			
                          dt=NULL,	    	# Sample time step
                          units=c("days", "hours", "mins", "secs"),		# Time units of sample time step
                          Nfilt=NULL,		  # Moving average filter width for detrending time series
                          nf=1,		# Number of moving average filter passes
                          NLowPass = NULL,	# Moving average filter width for low pass filtering demodulated signal
                          fDL=NULL,		    # Relative fraction of light hours during the day
                          MAcorrect = FALSE,	# Correction factor for moving average filter
                          phi= 51.176,
                          lambda=4.326,
                          ...)			      # Other parameters to be passed to O22GPpreprocess
{
  ts <- x # internal alias
  
  if (inherits(ts,"data.frame")){
    if (any(is.na(ts[,c(1,2)]))) stop("NAs not allowed")
    if (dim(ts)[2]>=2){
      t <- ts[,1]
      if (inherits(t, "POSIXt")){ 
        if (any(0!=diff(diff(round.POSIXt(t, units="secs"))))) stop("Samples are not equidistant and/or gaps are present in the series. Re-sample, fill gaps and/or split series in continuous subsets")    
        dt <- difftime(t[2],t[1],units="days")        
      } else {
        if (any(0!=round(diff(diff(t)), digits=5))) stop("Samples are not equidistant and/or gaps are present in the series. Re-sample, fill gaps and/or split series in continuous subsets")      
        dt <- diff(c(t[1],t[2]))  
      }
      if (is.null(fDL)){
        if (is.null(phi)|is.null(lambda)) {stop("Provide longitude and lattitude")} else {
        fDL <- fDLfun(t, phi=phi, lambda=lambda)
        fDLfilt <- filter(fDL, rep(1/Nfilt, Nfilt))
        }
      } else if (length(fDL)==1){
        fDL <- rep(fDL, length(ts[,2]))
      } else if (length(fDL)!= length(ts[,2])){
        stop("When fDL is provided, it should either be single valued or a vector of same length as the time series")
      }
      O2 <- ts[,2]
    }
  } else {
    if (is.null(dt)) stop("Provide dt")
    if (is.null(fDL)) stop("Provide daylight fraction of day (fDL)")
    if (any(is.na(ts))) stop("NAs not allowed")
    O2 <- ts
    fDL <- rep(fDL, length(O2))
    unit <- match.arg(units)
    if (!inherits(dt,"difftime")){
      if (is.null(units)) {
        stop("Either provide dt as difftime object, or provide units of dt")
      } else {
        dt <- as.difftime(dt, units=unit)
      }
    }
    t <- 1:length(O2)*dt
  }
  if (inherits(dt, "difftime")){
    dtnumeric <- as.double(dt, units="days")
  } else {
    dtnumeric <- dt
  }
  
  
  if (is.null(Nfilt)) stop("Provide moving average filter width (Nfilt)")
  if (Nfilt>length(O2)) stop("Filter longer than time series")
  
  # angular frequency (cycles per day) associated with diurnal periodicity, i.e. T = 1 day
  w <- 2*pi/1
  
  xpreprocessed <- GPPFourierPreprocess(O2, dt=dt, Nfilt=Nfilt,...)
  xfilt <- xpreprocessed$filt
  xres <- xpreprocessed$res
  
  timessub <- t[xpreprocessed$indices]
  fDLsub <- fDL[xpreprocessed$indices]
  
  # demodulate time series
  xdemod <- demod(xres,dt=dtnumeric,Nf=NLowPass, nf=nf)
  
  # correct for truncated sinusoid form of GPP(t)
  theta1 <- pi*(1-fDLsub)
  truncFac <- (sin(theta1) + (pi-theta1) * cos(theta1))/(pi-theta1 + 0.5 * sin(2*theta1))
  
  if (MAcorrect){
    # correct for moving average frequency response
    L <- Nfilt
    maxfreq <- (length(xdemod) - 1)/(length(xdemod))/dt
    w24h <- 2*pi/maxfreq
    MAcorrection <- 1/(1-Hw(w24h,L))
  } else{
    MAcorrection <- 1
  }
  
  # collect all correction factors
  xdemodcorrect <- xdemod*truncFac*Mod(MAcorrection)
  
  # time derivative in the frequency domain
  GPt <- w*Mod(xdemodcorrect)
  
  return(data.frame(times=timessub,GPPt=GPt))
}
