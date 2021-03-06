% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dispersion.R
\name{select_variables}
\alias{select_variables}
\alias{select_rows}
\title{Select variables based on mean, var, and dispersion}
\usage{
select_variables(
  dat,
  meta = NULL,
  filter.zdisp = NULL,
  filter.disp = NULL,
  filter.var = NULL,
  filter.mean = NULL,
  n_bins = 10,
  verbose = FALSE,
  seed = NULL,
  ...
)

select_rows(
  dat,
  meta = NULL,
  filter.zdisp = NULL,
  filter.disp = NULL,
  filter.var = NULL,
  filter.mean = NULL,
  n_bins = 10,
  verbose = FALSE,
  seed = NULL,
  ...
)
}
\arguments{
\item{dat}{a time-series data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).}

\item{meta}{a data frame consisted of meta information about m variables. If given, this data frame will be thresholded accordingly.}

\item{filter.zdisp}{a range of z-scores of dispersion to select the variables. Variables with z-scores of dispersion that are outside of this range will be excluded.}

\item{filter.disp}{a range of dispersion values to select the variables. Variables with dispersions that are outside of this range will be excluded.}

\item{filter.var}{a range of variances to select the variables. Variables with variances that are outside of this range will be excluded.}

\item{filter.mean}{a range of means to select the variables. Variables with means that are outside of this range will be excluded.}

\item{verbose}{a logical specifying to print the computational progress. By default, \code{FALSE}.}

\item{seed}{a seed for the random number generator.}

\item{...}{optional arguments.}

\item{data.only}{a logical specifying to return only the new data matrix with selected variables. By default, \code{TRUE}.}
}
\value{
If \code{data.only=TRUE}, \code{select_variables} returns a new data matrix containing only the selected variables.

If \code{data.only=FALSE}, \code{select_variables} returns a list, consisted of a new data matrix and a data frame of calculated statistics for variables (see \code{dispersion}).
}
\description{
This function calculates the mean, variance, dispersion, and z-score of normalized dispersion for each variable (row).
Then, select variables that meet all of criteria given by filter.disp, filter.zdisp, filter.var, and filter.mean.
It's highly recommended to investigate the distribution of these statistics, because making thresholds.
}
\examples{
\dontrun{
data(cys_optm_missing)
meta <- cys_optm_missing[,1:3]
optm <- log(cys_optm_missing[,4:9])
days <- as.numeric(colnames(optm_missing))

disp_optm <- dispersion(optm)
disp_optm <- cbind(meta, disp_optm)
# make a histogram of dispersion statistics
hist(disp_optm$disp, 100)
# make a histogram of z-score of normalized dispersion
hist(disp_optm$zdisp, 100)

select_variable(optm, filter.zdisp = c(-2,2))
# library(readr)
# write_excel_csv(disp_optm, file="~/coptm_dispersion.csv")
}
}
\author{
Neo Christopher Chung \email{nchchung@gmail.com}
}
