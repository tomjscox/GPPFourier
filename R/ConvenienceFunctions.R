#
# 2014-12-12 Time Series Convenience functions
# 


#' Subsampling, interpolating and gapfilling
#' 
#' 



#' @title Interpolate, subsample and fill gaps in time series by linear interpolation
#' 
#' @description Interpolate, subsample and fill gaps in time series by linear interpolation
#' 
#' 
#' @param x Data frame with time in first column and time series data in
#' second. 
#' @param dt Sampling interval. Default to time difference of first two observations
#' @param units Optional: units of \code{dt}. Only necessary when time is of class \code{POSIXt} and \code{dt} is not of class \code{difftime}. 
#' @param type Type of interpolation. Only linear interpolation implemented
#' @return A data frame with time and interpolated time series.
#' @description Subsample and interpolate oversampled or undersampled time series to obtain equidistant data at a lower or higher sampling frequency. Gaps are filled with linear interpolation
#' @details An equidistant times series with sampling interval \code{dt} is created by linear interpolation. \code{gapfill} takes the time interval between the first and second observation as default value for \code{dt}. \code{subsample} and \code{resample} require a value of \code{dt}. 
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @keywords utilities
#' 
#' @export gapfill
gapfill <- function(x, 
                    dt=x[2,1] - x[1,1], 
                    units=c("auto", "secs", "mins", "hours","days", "weeks"), 
                    type="linear") {
  
  if (!type=="linear") stop("Only linear interpolation implemented")
  if (inherits(x[,1], "POSIXt")&!inherits(dt, "difftime")) 
  {
    if (is.null(units)){
      stop("Either dt must be of class difftime or units must be provided")  
    } else {
      dt <- as.difftime(dt,units=match.arg(units))
    }
  }
  
  
  t <- x[,1]
  y <- x[,2]
  
  if (length(y)>1){
    # dt <- t[2] - t[1]
    times <- seq(t[1],t[length(t)],dt)
    
    if(all(is.na(y))){
      yfilled <- rep(NA, length(times))
    } else {
      yfilled <- approx(t,y,times)$y 
    }
    
  } else {
    yfilled <- y
    times <- t
    warning("Only one observation, nothing done")
  }
  return(data.frame(t=times,y=yfilled))
}


#' @rdname gapfill
#' @export subsample
subsample <- function(x, dt, units=c("auto", "secs", "mins", "hours","days", "weeks")){
  # subsample(x=x,dt=dt, units=units)
  if (is.null(dt)) stop("dt has to be provided")
  gapfill(x=x, dt=dt, units=units)
}


#' @rdname gapfill
#' @examples 
#' plot(Hoernum, type="p", xlim=as.POSIXct(c("2008-08-01","2008-09-30")))
#' points(gapfill(Hoernum), type="p", pch=20, col="red", cex=0.2)
#' @export resample
resample <- function(x, dt, units=c("auto", "secs", "mins", "hours","days", "weeks")){
# subsample(x=x,dt=dt, units=units)
  if (is.null(dt)) stop("dt has to be provided")
  gapfill(x=x, dt=dt, units=units)
}




# 2013-05-22 add functionality to return the trend. Should be backward compatible
# Old version of detrend below


#' @title Detrend a vector
#' 
#' @description Substract a linear trend from a vector
#' 
#' @details The linear trend is calculated as \code{lm(x~1:length(x))}
#' 
#' @param x Vector or time series to be detrended
#' @param returntrend When FALSE the detrended vector is
#' returned, when TRUE a 2 element list is returned with both the trend and the
#' detrended vector
#' @return The detrended vector. When \code{returntrend = TRUE} a 2 element
#' list is returned with the trend and the detrended vector.
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @keywords utilities
#' @export detrend
detrend <- function(x, returntrend = FALSE){
t <- 1:(length(x))
trendmod <- lm(x~t)
trend <- predict(trendmod)

if (returntrend) {
return(list(detrended=x-trend, trend=trend))
} 
else {
return(x-trend)
}
}

# Old version of detrend
# detrend <- function(y, returntrend = FALSE){
# x <- 1:(length(y))
# trendmod <- lm(y~x)
# trend <- predict(trendmod)
# return(y-trend)
# }
