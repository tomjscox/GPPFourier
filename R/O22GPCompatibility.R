# 2016-09-09 Manuscript was written with older version of GPPFourier package. Functions have since been renamed from O22GPxxxx to GPPFourierxxxx.
# GPPFourier package compatibility functions
# With these functions no changes in Manuscript processing scripts is necessary to use GPPFourier package in stead of the old O22GP package

#' @rdname O22GP
#' @export O22GPpreprocess
O22GPpreprocess <- function(...){
args <- list(...)
removeargs <- match(c("bw", "bwfreq"), names(args))
removeargs <- removeargs[!is.na(removeargs)]
args <- args[-removeargs]
# args
do.call(GPPFourierPreprocess,args)
}



#' @title Deprecated: Calculate GPP from O2 time series
#' 
#' @description Don't use this function, use GPPFourier() instead. This function only exists
#' for backward compatibility of GPPFourier package with older scripts.
#' 
#' 
#' @aliases O22GP WindowO22GP WindowO22GP.gts
#' @param ... All arguments are passed to equivalent GPPFourier functions
#' @return
#' 
#' See \code{\link{GPPFourier}}, \code{\link{GPPFourierPreprocess}},
#' \code{\link{WindowGPPFourier}}, \code{\link{WindowGPPFourier.gts}}
#' @author Tom Cox <tom.cox@uantwerp.be>
#' @keywords utilities
#' @export O22GP
O22GP <- function(...){
GPPFourier(...)
}

#' @rdname O22GP
#' @export WindowO22GP
WindowO22GP <- function(...){
WindowGPPFourier(...)
}

#' @rdname O22GP
#' @export WindowO22GP.gts
WindowO22GP.gts <- function(...){
WindowGPPFourier.gts(...)
}

