#' Calculate dispersion metrics
#'
#' For each variable (row), calculate a dispersion (var/mean) and a z-score of normalized dispersion.
#' A z-score of normalized dispersion is motivated due to correlation between means and variances.
#' Variables are binned (by default,\code{n_bins=10}) by their means. Within a bin, a z-score is calculated.
#' It's recommended that you look at both results \code{zdisp} and \code{disp} to decide how to filter the data.
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
#' @return \code{dispersion} returns a data.frame of statistics for \code{m} variables (rows).
#'
#' @export dispersion
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
#' disp_optm <- dispersion(optm, timepoints = days, dof="cv")
#' disp_optm <- cbind(meta, disp_optm)
#' # make a histogram of dispersion statistics
#' hist(disp_optm$disp, 100)
#' # make a histogram of z-score of normalized dispersion
#' hist(disp_optm$zdisp, 100)
#' # library(readr)
#' # write_excel_csv(disp_optm, file="~/coptm_dispersion.csv")
#'}
dispersion <- function(dat,
                       n_bins = 10,
                        center.dat = TRUE,
                        scale.dat = FALSE,
                        verbose = FALSE,
                        seed = NULL,
                        ...) {

  if (is.null(seed)) set.seed(seed)
  m <- nrow(dat)
  n <- ncol(dat)
  if(scale.dat | center.dat) dat <- t(scale(t(dat),center=center.dat,scale=scale.dat))

  rmean <- apply(dat, 1, function(x) mean(x,na.rm = TRUE))
  rvar <- apply(dat, 1, function(x) var(x,na.rm = TRUE))
  rdisp <- rvar/rmean
  if(verbose) plot(log(rmean),log(rvar),main="Relationship between Means and Variances",pch=20)

  rcut <- cut(rmean, breaks=n_bins)
  rzdisp <- vector("numeric", length=m)
  for(ircut in unique(rcut)) {
    if(sum(rcut == ircut) == 1) {
      rzdisp[rcut == ircut] = rdisp[rcut == ircut]
    } else {
      rzdisp[rcut == ircut] <- (rdisp[rcut == ircut] - mean(rdisp[rcut == ircut])) / sd(rdisp[rcut == ircut])
    }
  }

  #### RETURN
  output <- data.frame(mean=rmean, var=rvar, disp=rdisp, zdisp=rzdisp)
  return(output)
}

#' Select variables based on mean, var, and dispersion
#'
#' This function calculates the mean, variance, dispersion, and z-score of normalized dispersion for each variable (row).
#' Then, select variables that meet all of criteria given by filter.disp, filter.zdisp, filter.var, and filter.mean.
#' It's highly recommended to investigate the distribution of these statistics, because making thresholds.
#'
#' @param dat a time-series data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param meta a data frame consisted of meta information about m variables. If given, this data frame will be thresholded accordingly.
#' @param filter.disp a range of dispersion values to select the variables. Variables with dispersions that are outside of this range will be excluded.
#' @param filter.zdisp a range of z-scores of dispersion to select the variables. Variables with z-scores of dispersion that are outside of this range will be excluded.
#' @param filter.mean a range of means to select the variables. Variables with means that are outside of this range will be excluded.
#' @param filter.var a range of variances to select the variables. Variables with variances that are outside of this range will be excluded.
#' @param data.only a logical specifying to return only the new data matrix with selected variables. By default, \code{TRUE}.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param ... optional arguments.
#'
#' @return If \code{data.only=TRUE}, \code{select_variables} returns a new data matrix containing only the selected variables.
#' @return If \code{data.only=FALSE}, \code{select_variables} returns a list, consisted of a new data matrix and a data frame of calculated statistics for variables (see \code{dispersion}).
#'
#' @export select_variables
#' @import splines
#' @import stats
#' @import graphics
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'@examples
#'\dontrun{
#' data(cys_optm_missing)
#' meta <- cys_optm_missing[,1:3]
#' optm <- log(cys_optm_missing[,4:9])
#' days <- as.numeric(colnames(optm_missing))
#'
#' disp_optm <- dispersion(optm)
#' disp_optm <- cbind(meta, disp_optm)
#' # make a histogram of dispersion statistics
#' hist(disp_optm$disp, 100)
#' # make a histogram of z-score of normalized dispersion
#' hist(disp_optm$zdisp, 100)
#'
#' select_variable(optm, filter.zdisp = c(-2,2))
#' # library(readr)
#' # write_excel_csv(disp_optm, file="~/coptm_dispersion.csv")
#'}
select_variables <- function(dat,
                       meta = NULL,
                       filter.zdisp = NULL,
                       filter.disp = NULL,
                       filter.var = NULL,
                       filter.mean = NULL,
                       n_bins = 10,
                       verbose = FALSE,
                       seed = NULL,
                       ...) {

  if (is.null(seed)) set.seed(seed)
  m <- nrow(dat)
  n <- ncol(dat)
  if(is.null(filter.zdisp) & is.null(filter.disp) & is.null(filter.var) & is.null(filter.mean)) stop("None of the filters are specified.")

  dat_disp <- dispersion(dat)

  select <- apply(dat_disp, 1, function(x) {
        if(!is.null(filter.zdisp)) tmp <- (x[4] > filter.zdisp[1]) & (x[4] < filter.zdisp[2]);
        if(!is.null(filter.disp)) tmp <- tmp & (x[3] > filter.disp[1]) & (x[3] < filter.disp[2]);
        if(!is.null(filter.var)) tmp <- tmp & (x[2] > filter.var[1]) & (x[2] < filter.var[2]);
        if(!is.null(filter.mean)) tmp <- tmp & (x[1] > filter.mean[1]) & (x[1] < filter.mean[2]);
        return(tmp)
      })

  #### RETURN
  if(data.only) {
    return(dat[select,])
  } else {
    dat.selected <- dat[select,]
    if(!is.null(meta)) meta <- meta[select,]
    output <- list(dat.selected = dat, dat.original=dat, meta = meta, dispersion = dat_disp)
    return(output)
  }
}

#' @rdname select_variables
#' @export
select_rows <- select_variables
