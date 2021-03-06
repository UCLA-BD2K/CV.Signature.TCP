#' Preprocess the input data by fitting cubic splines
#'
#' The data is denoised and/or imputed by fitting cubic splines. This requires a minimum of 4 data points per variables (e.g., proteins).
#' If a degree of freedom (\code{dof}) is set to \code{cv}, cross validation is performed on each variable to identify the optimal dof for that variable.
#' If \code{dof} set to \code{cv.global}, the mean of all cross-validated DoFs is used for all variables.
#' Lastly, \code{dof} can be set to a fixed numeric value predetermied by the user.
#'
#' @param dat a time-series data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param timepoints a vector of time points for columns of dat.
#' @param dof the degree of freedom in spline.smooth. By default, \code{cv} performs cross-validation per-variable.
#' @param center.dat a logical specifying to center the input and denoised data. By default, \code{TRUE}.
#' @param scale.dat a logical specifying to scale the input and denoised data. By default, \code{FALSE}.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param ... optional arguments.
#'
#' @return \code{preprocess_spline} returns an imputed data, a denoised and imputed data data, and a vector of degrees of freedom.
#'
#' @export preprocess_spline
#' @import splines
#' @import stats
#' @import graphics
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'@examples
#'\dontrun{
#' data(cys_optm)
#' meta <- cys_optm[,1:4]
#' optm <- log(cys_optm[meta$Select,5:10])
#' optm <- t(scale(t(optm), scale=TRUE, center=TRUE))
#' days <- as.numeric(colnames(optm))
#'
#' preprocessed_optm <- preprocess_spline(optm, timepoints = days, dof="cv")
#'}
preprocess_spline <- function(dat,
                          timepoints=NULL,
                          dof=c("cv","cv.global"),
                          center.dat = TRUE,
                          scale.dat = FALSE,
                          verbose = FALSE,
                          seed = NULL,
                          ...) {

  if (is.null(seed)) set.seed(seed)
  m <- nrow(dat)
  n <- ncol(dat)
  if(scale.dat | center.dat) dat <- t(scale(t(dat),center=center.dat,scale=scale.dat))

  # sanity check
  if(n != length(timepoints)) stop("denoise_spline: The number of time points must match the number of columns in the input data.")

  dat.denoise = matrix(0, nrow=m, ncol=n)
  df <- vector("numeric",m)
  if(dof == "cv") {
    if(verbose) message("\n\r Individual degrees of freedom via cross validation:")
    for(i in 1:m) {
      if(verbose) cat(i)

      #fitting smoothing splines using smooth.spline(X,Y,df=...)
      tpi <- data.frame(x=timepoints, y=dat[i,])
      spline.fit <- with(tpi[!is.na(tpi$y),], smooth.spline(x,y,cv = TRUE))
      df[i] <- spline.fit$df
      dat.denoise[i,] <- with(tpi, predict(spline.fit, x)$y)
    }
  } else if(dof == "cv.global") {
    if(verbose) message("\n\r Cross Validation:")
    for(i in 1:m) {
      if(verbose) cat(i)
      tpi <- data.frame(x=timepoints, y=dat[i,])
      df[i] <- with(tpi[!is.na(tpi$y),], smooth.spline(x,y,cv = TRUE)$df)
    }
    if(verbose)  message("\n\r Global degree of freedom:")
    for(i in 1:m) {
      if(verbose) cat(i)
      tpi <- data.frame(x=timepoints, y=dat[i,])
      globaldf = median(df, na.rm = TRUE)
      spline.fit <- with(tpi[!is.na(tpi$y),], smooth.spline(x,y,df=globaldf))
      dat.denoise[i,] <- with(tpi, predict(spline.fit, x)$y)
    }
    if(verbose) {
      message(paste("The median DoF is used for all variables:",globaldf))
      hist(df,20,col="grey",main="Degrees of Freedom from Cross Validation")
      abline(v=globaldf,col="red",lty="dashed")
    }
  } else if(is.numeric(dof)){
    df <- dof
    if(verbose) message("\n\r The degree of freedom (dof) is set to the user input:")
    for(i in 1:m) {
      if(verbose) cat(i)

      #fitting smoothing splines using smooth.spline(X,Y,df=...)
      spline.fit <- smooth.spline(x=timepoints,y=dat[i,],df=dof)
      dat.denoise[i,] <- predict(spline.fit,newdata = timepoints)$y
    }
  }
  rownames(dat.denoise) <- rownames(dat)
  colnames(dat.denoise) <- colnames(dat)

  dat.impute <- dat
  dat.impute[is.na(dat.impute)]  <- dat.denoise[is.na(dat.impute)]

	return(list(dat.impute=dat.impute,
				dat.denoise=dat.denoise,
				imputed=is.na(dat),
				df=df))
}

#' @rdname denoise_pca
#' @export
denoise_spline <- preprocess_spline
