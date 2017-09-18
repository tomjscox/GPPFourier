#
# Test function with enhanced functionality
#


#' GPPFourier test function
#' 
#' GPPFourier test function
#' 
#' This function allows to play around with the basics behind GPPFourier and
#' GPPFourierPreprocess
#' 
#' \code{trunclight} controls whether GPP is assumed to evolving as a truncated
#' sinusoid over a day or as a pure sinusoidal.  This assumption determines the
#' relation between time averaged GPP and the Fourier amplitude at diel
#' frequency
#' 
#' \code{usefft}: by default \code{GPPFourier()} calculates the amplitude at
#' diel frequency directly \code{fft()} calculates the full Fourier transform
#' of the time series. The amplitude at diel frequency can be derived from the
#' Fourier transform
#' 
#' \code{fourierderive}: if TRUE time derivation is performed in the Fourier
#' domain by multiplying with i x omega, omega being the diurnal angular
#' frequency.  Thus the GPPFouriertest returns the amplitude at diel frequency of dx/dt. If
#' fourierderive is FALSE, the amplitude at diel frequency of x is returned
#' (multiplied with the factor determined by \code{trunclight}
#' 
#' @param x O2 time series from which GPP is to be calculated
#' @param dt Sample time step
#' @param units Time unit of sampling time step
#' @param Detrend Toggle time series detrending. See GPPFourierPreprocess
#' @param filter Toggle time series filtering. See GPPFourierPreprocess
#' @param Nfilt Moving average filter width
#' @param circular Moving average circular boundary condition
#' (see documentation of filter())
#' @param sides Moving average central or one sided (see documentation of
#' filter())
#' @param filtcorrect Logical controlling whether GPP estimate
#' is corrected for signal falsely removed by filtering
#' @param trunclight Use truncated sinusoid approximation for light&GPP. If
#' FALSE, a sinusoid approximation is assumed
#' @param fDL Relative fraction of light hours during the day
#' @param usefft If FALSE the amplitude at diel freqyency is computed directly. If TRUE fft() is used to estimate amplitude at diel frequency. 
#' @param padlength Number of zeroes to be appended to the time series to
#' increase frequency resolution
#' @param fourierderive Calculate derivative in the frequency domain or not
#' @param confine Confine time series to integer number of days or tidal cycles
#' @param taper Taper the time series with spec.taper()
#' @param p Fraction of the time series to be tapered at each side
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @references Cox T.J.S. et al. (2015) 'Estimating primary production from oxygen time series: a novel approach in the frequency domain', Limnology And Oceanography:Methods 13, 529-552. DOI: 10.1002/lom3.10046

#' @seealso \link{GPPFourier}, \link{GPPFourierPreprocess}
#' @keywords utilities
#' @export GPPFouriertest
GPPFouriertest <- function(	x,		# O2 time series from which GPP is to be calculated
			dt=NULL,		# Sample time step
			units="days",		# Time units of sample time step
			Detrend=FALSE,		# Toggle time series detrending
			filter = FALSE,		# Toggle time series filtering
			Nfilt=NULL,		# Moving average filter width
			circular = FALSE,	# Moving average circular boundary condition (see documentation of filter())
			sides = 2,		# Moving average central or one sided (see documentation of filter())
			filtcorrect = TRUE,	# Correct for frequency response of moving average filter 
			trunclight = TRUE,	# Use truncated sinusoid approximation for light&GPP. If FALSE, a sinusoid approximation is assumed
			fDL=NULL,		# Relative fraction of light hours during the day
			usefft = FALSE,		# Use fft to estimate peak corresponding with diurnal periodicity.  
			padlength=0,		# Number of zeroes to be appended to the time series to increase frequency resolution
			fourierderive = FALSE,	# Calculate derivative in the frequency domain or not
			confine = "days",	# Confine time series to integer number of days or tidal cycles 
			taper = FALSE,		# Taper the time series with spec.taper()
			p = 0.1)		# Fraction of the time series to be tapered at each side
			
{
# Test function with enhanced functionality. 
# Calculate average gross production from O2 time series
# Units: unit of oxgen as in ts per time unit (default = days)
#

if (is.null(dt)) stop("Provide dt")
if (filter&is.null(Nfilt)) stop("Provide filter length (Nfilt)")
if (trunclight&is.null(fDL)) stop("Provide daylight fraction of day (fDL)")
if (Nfilt>length(x)) stop("Filter longer than time series")

#Frequency response of moving average filter
Hw <- function(w,L){
1/L*(1-exp(-1i*L*w))/(1-exp(-1i*w))
}

# angular frequency (cycles per day) associated with diurnal periodicity, i.e. T = 1 day
w <- 2*pi/1

# Pre-process time series
xPreProcess <- GPPFourierPreprocess(x, dt=dt, units=units, Detrend=Detrend, filter = filter, Nfilt=Nfilt, circular = circular, sides = sides)

filtx <- xPreProcess$filt
resx <- xPreProcess$res
trendx <- xPreProcess$trend

if (taper){
resx <- spec.taper(resx,p=p)
}

if (usefft) {
  fx <- fft(c(resx,rep(0,padlength)))/length(x)
  dfres <- 1/(length(resx)+padlength)/dt
  freq <- 0:(length(resx)+padlength-1)*dfres
  i24h <- round(1/df+1)
  A24 <- 2*fx[i24h]
} else {
  times <- (1:length(resx))*dt
  A24 <- 2*sum(exp(-times*1i*w)*resx)/length(resx)
}

if (trunclight){
theta1 <- pi*(1-fDL)
truncFac <- (sin(theta1) + (pi-theta1) * cos(theta1))/(pi-theta1 + 0.5 * sin(2*theta1))
A24 <- A24*truncFac
} else {truncFac <- 1}

if (fourierderive){
A24 <- 1i*w*A24
}

if (filtcorrect) {
# Properties moving average filter filter
L <- Nfilt
maxfreq <- (length(resx) + padlength - 1)/(length(resx) + padlength)/dt
w24h <- 2*pi/maxfreq
A24 <- A24/(1-Hw(w24h,Nfilt))
#omega <- (0:3000)/100Fx[i24h]0
#plot(omega,Mod(Hw(omega,L)),type="l",xlim=c(0,3))
#plot(omega,Arg(Hw(omega,L)),type="l",xlim=c(0,.1))
#Mod(Hw(w24h,Nfilt))
#Arg(Hw(w24h,Nfilt))
}

return(A24)

}
