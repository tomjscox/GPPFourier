# Auxiliary functions used in processing O2-data
# 2013-05-22 extended functionality of detrend(). 
# 2013-07-30 added plotF from project artikel-fourier/R/functions/FourierFunctions.R
# 2014-12-12 added High level functions WindowGPPFourier.gts and WindowGPPFourier


#
# Calculate fft, padded with zeroes and plot result.
#

#
# Increase resolution of fft by padding zeroes behind the series
#



#' @title Fourier transform of padded time series
#' 
#' @description Fourier transform of padded time series
#' 
#' 
#' @param x Time series for which fourier transform (fft) is calculated after
#' padding with zeroes
#' @param n Number of zeroes to be padded
#' @return Fourier transform of x padded with n zeroes
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @keywords utilities
#' # @export paddedfft
#' @keywords internal
paddedfft<- function(x, # series to be fourier transformed (fft)
                     n # number of zeros to be padded
)
{
  x2 <- x - mean(x)
  x2<- c(x2,rep(0,n))
  fx2 <- fft(x2)
  return(fx2 # Fourier transform of x padded with n zeroes
  )
}




#' Test plot for padded fft calculation
#' 
#' Test plot for padded fft calculation
#' 
#' 
#' @param y Vector or time series for which fft and padded fft is calculated
#' @param padlength Number of zeroes to be appended to \code{y}
#' @param dt Time interval between succesive values of \code{y}
#' @param ... Other arguments passed to plot()
#' @return A two panel plot is generated with the unpadded fft on top and the
#' padded fft at the bottom
#' @author Tom Cox <tom.cox@@uantwerp.be>
#' @keywords utilities
#' #@export padTestPlot
#' @keywords internal
padTestPlot <- function(y,padlength,dt=1,...) {

fy <- fft(y-mean(y))
df = 1/dt
f <- ((1:length(y))-1)/length(y)*df

padfy <- fft(c(y-mean(y),rep(0,padlength)))
paddf = 1/(length(y)+padlength)/dt
padf <- ((1:length(padfy))-1)*paddf

par(mfcol=c(2,1))
plot(f,Mod(fy)/length(y),xlim=c(0,5),type="b",main="fft",...)
plot(padf,Mod(padfy)/length(y),xlim=c(0,5),type="b",main=paste("padded fft,n=",padlength),...)

}




#' Plot amplitude or argument of Fourier transform of a time series
#' 
#' Plot amplitude or argument of Fourier transform of a time series
#' 
#' 
#' @param x Time series
#' @param dt Time time interval, used to calculate frequencies
#' @param units Time unit, default = "days"
#' @param xlab x-label of plot. Default: f [d-1]" "f [cycles per day]"
#' @param ylab y-label of plot. Default: |F(x)|, or Arg(F(x)) when argument=TRUE
#' @param detrend If TRUE, series is detrended first. Default=FALSE
#' @param add If TRUE, plot is added to existing plot
#' @param argument If TRUE, plot the argument in stead of the modulus
#' @param ... Arguments passed to plot()
#' @return A plot of the amplitude of the Fourier transform, calculated with
#' fft(). No zero-padding or tapering is performed.
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @keywords utilities
#' @importFrom  graphics lines par plot points
#' @importFrom stats approx df fft filter lm predict spec.taper
#' @export plotF
plotF <- function(x, 					#Time series 
			dt = 1,				#Sampling time interval, used to calculated frequencies
			units = "days",
			xlab=NULL,
			ylab=NULL,			# "|F(x)|"
			detrend=FALSE,			#
			add = FALSE,
			argument = FALSE,	 	# If TRUE, plot the argument in stead of the modulus
			...)				# Arguments passed to plot

{
# Plots the Fourier transform with meaningfull frequencies

if (detrend==TRUE){
N <- length(x)
index <- 1:N
x <- x - predict(lm(x~index))						#Detrend series		
}

if (is.null(ylab)){
if (argument){
  ylab <- "Arg(F(x))"
} else {
  ylab <- "|F(x)|"
}
}
  
if (is.null(xlab)){
  xlab <- switch(units,
         "days"=expression(paste("f [d"^-1,"]")),
         "mins"=expression(paste("f [min"^-1,"]")),
         "secs"=expression(paste("f [s"^-1,"]")))
}

Fx <- fft(x)
df = 1/dt
f <- ((1:length(x))-1)/length(x)*df
if (add==FALSE){
if (!argument){
plot(f,Mod(Fx)/length(x),xlab=xlab,ylab=ylab,...)
} else{
plot(f,Arg(Fx),xlab=xlab,ylab=ylab,...)
points(f,Arg(Fx)+2*pi,xlab=xlab,ylab=ylab,...)
points(f,Arg(Fx)-2*pi,xlab=xlab,ylab=ylab,...)
}
} else {
if (!argument){
lines(f,Mod(Fx)/length(x),...)
} else {
lines(f,Arg(Fx),...)
points(f,Arg(Fx)+2*pi,xlab=xlab,ylab=ylab,...)
points(f,Arg(Fx)-2*pi,xlab=xlab,ylab=ylab,...)
}
}
return(invisible(data.frame(w=f,F=Fx/length(x))))
}



#' Frequency response of moving average filter
#' 
#' Frequency response of moving average filter
#' 
#' 
#' @param w The frequency
#' @param L Length of the moving average filter
#' @return Hw(w,L) = 1/L*(1-exp(-1i*L*w))/(1-exp(-1i*w))
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @keywords utilities
#' @export Hw
Hw <- function(w,L){
1/L*(1-exp(-1i*L*w))/(1-exp(-1i*w))
}
