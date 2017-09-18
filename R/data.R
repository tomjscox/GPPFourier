
#
# 2016-09-13 This file was added in the R-directory with only purpose of generating Rd files with roxygen2 
#


#' @title Example time series 
#' @description Time series observed in 2010 in the Scheldt estuary (Kruibeke pontoon, N51.176 E4.326) and at Hoernum Tief in 2008 in the German Wadden Sea (Hoernum Tief pole, N54.783 E8.45). 
#' Depth averaged O2 and GPP time series, simulated with a 1D vertical model of O2 dynamics. O2 and GPP time series at a fixed estuarine location simulated with 1D longitudinal estuary model of O2 dynamics  (Cox et al, 2017)
#' @name Kruibeke
# #' @aliases Hoernum, watercolumn, estuary
#' @docType data
#' @format A data frame containing observation time, the observed dissolved oxygen concentration in \if{latex}{\out{$\mu$}} \if{html}{\out{&mu;}}M
#' @references 
#' Cox T.J.S. et al. (2015) 'Estimating primary production from oxygen time series: a novel approach in the frequency domain', Limnology And Oceanography:Methods 13, 529-552. DOI: 10.1002/lom3.10046
#' Cox T.J.S. et al. (2017) 'Tune in on 11.57 \if{latex}{\out{$\mu$}} \if{html}{\out{&mu;}}Hz and listen to primary production', Biogeosciences Discussions doi:10.5194/bg-2017-81
#' @keywords datasets
#' @examples
#' plot(Kruibeke, pch=20)
#' plot(Hoernum, pch=20)
#' plot(watercolumn[,c("time", "O2")], pch=20)
#' plot(estuary[,c("time", "O2")], pch=20)
NULL

#' @name watercolumn
#' @rdname Kruibeke
NULL

#' @name estuary
#' @rdname Kruibeke
NULL

#' @name Hoernum
#' @rdname Kruibeke
NULL

