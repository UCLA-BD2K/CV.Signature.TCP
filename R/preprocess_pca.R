#' Preprocess the input data by principal component analysis (PCA) and related methods
#'
#' The data is denoised and the missing values are imputed by using the top r eigenmatrices.
#' After apply SVD on the centered and/or scaled data, the top r eigenmatrices are constructed where is r < min(n,m).
#'
#' To impute missing values with PCA/SVD, two approximation methods are provided.
#' For method="nipals", a Non-linear Iterative Partial Least Squares (NIPALS) algorithm is used from \code{nipals} in the mixOmics package.
#' For method="em", a low-rank SVD approximation by the EM algorithm is used from \code{imputed.svd} in the bcv package.
#'
#' @param dat a time-series data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param r the number of PCs or eigenmatrices to retain.
#' @param method a method to perform singular-value decomposition when a dataset has missing values. See below for the explanation.
#' @param center.dat a logical specifying to center the input and denoised data. By default, \code{TRUE}.
#' @param scale.dat a logical specifying to scale the input and denoised data. By default, \code{FALSE}.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param ... optional arguments.
#'
#' @return \code{preprocess_pca} returns a matrix of imputed and/or denoised data.
#'
#' @export preprocess_pca
#' @importFrom mixOmics nipals
#' @importFrom corpcor fast.svd
#' @importFrom bcv impute.svd
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T., Tibshirani, R., Botstein, D. and Altman, R.B. (2001). Missing value estimation methods for DNA microarrays. Bioinformatics 17(6), 520–525.
#'@examples
#'\dontrun{
#' data(cys_optm)
#' meta <- cys_optm[,1:4]
#' optm <- log(cys_optm[meta$Select,5:10])
#' optm <- t(scale(t(optm), scale=TRUE, center=TRUE))
#' days <- as.numeric(colnames(optm))
#'
#' preprocessed_optm <- preprocess_pca(optm, timepoints = days, r=3, method=c("nipals","em"))
#'}
preprocess_pca <- function(dat,
                r=NULL,
                method = c("em","nipals"),
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
  if(r >= ncol(dat) | r >= nrow(dat)) stop("preprocess_pca: The number of PCs r < min(n,m).")

  if(is.numeric(r)) {
    if(verbose) message(paste("The number of PCs to retain:",r))
  } else {
    stop("\n\r To denoise the data with PCA, r must be a numeric value corresponding to the number of PCs to retain.")
  }

  if(method=="em") {
    dat.denoise <- bcv::impute.svd(dat, k=r, ...)$x
  } else if(method=="nipals") {
    dat.denoise <- mixOmics::nipals(dat, reconst = TRUE, ncomp = r, ...)$rec
  }

  rownames(dat.denoise) <- rownames(dat)
  colnames(dat.denoise) <- colnames(dat)

  dat.impute <- dat
  dat.impute[is.na(dat.impute)]  <- dat.denoise[is.na(dat.impute)]

	return(list(dat.impute=dat.impute,
				dat.denoise=dat.denoise,
				imputed=is.na(dat)))
}

#' @rdname denoise_pca
#' @export
denoise_pca <- preprocess_pca
