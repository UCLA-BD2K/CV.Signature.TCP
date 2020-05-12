#' Wrapper to clustering algorithms
#'
#' This is a convenient wrapper to Kmeans and Hierarchical clustering when using the `tms` package.
#'
#' @param dat a (denoised) data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param timepoints a vector of time points for columns of dat. While this is not necessary, it's here for now to remind the user about the analysis pipeline.
#' @param K a number of clusters.
#' @param dist.method a distance method for time course data, resulting in a \code{m * m} distance matrix for rows. 'dtw' for dynamic time wrapping or 'cor.diss' for correlation-based dissimilarities.
#' @param cluster.method a clustering method.
#' @param hclust.algorithm an algorithm used in hierarchical clustering. See \code{hclust} for a full list.
#' @param kmeans.algorithm an algorithm used in kmeans clustering. See \code{kmeans} for a full list.
#' @param center.dat a logical specifying to center the input and denoised data. By default, \code{TRUE}.
#' @param scale.dat a logical specifying to scale the input and denoised data. By default, \code{FALSE}.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param ... optional arguments.
#'
#' @return \code{cluster} returns a list consisting of
#' \item{dat.dtw}{\code{m * m} distance matrix based on dynamic time wrapping.}
#' \item{cluster.obj}{an object returned from clustering the predicted data.}
#' \item{membership}{a vector of length \code{m}, identities of clusters.}
#'
#' @export cluster
#' @importFrom dtw dtwDist
#' @importFrom TSclust diss
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'@examples
#'\dontrun{
#' ## load the example "coptm" data
#' data(cys_optm)
#' meta <- cys_optm[,1:4]
#' optm <- log(cys_optm[meta$Select,5:10])
#' optm <- t(scale(t(optm), scale=TRUE, center=TRUE))
#' days <- as.numeric(colnames(optm))
#'
#' ## denoise using the cubic splines
#' denoised_optm <- denoise_spline(optm, timepoints = days, dof="cv", verbose=FALSE)
#'
#' ## cluster the denoised data using K-means clustering
#' clustered_optm <- cluster(denoised_optm,
#'                 timepoints = days,
#'                 cluster.method = "kmeans",
#'                 K=6,
#'                 center.dat = TRUE,
#'                 scale.dat = FALSE,
#'                 verbose = TRUE)
#'}
cluster <- function(dat,
                timepoints=NULL,
                K,
                dist.method = c("euclidean", "cor.diss", "dtw"),
                cluster.method = c("kmeans", "hclust"),
                hclust.algorithm = "complete", # see hclust for a full list
                kmeans.algorithm = "Hartigan-Wong", # see kmeans for a full list
                kmeans.centers,
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
  if(n != length(timepoints)) stop("The number of time points must match the number of columns in the input data.")
  if(cluster.method == "kmeans") {dist.method == "euclidean"}
  if(cluster.method=="kmeans" & dist.method=="dtw") { stop("You have chosen `dtw` and `kmeans`. To cluster a distance matrix from DTW, please use hclust.") }

  #### DISTANCE
  if(dist.method == "dtw") {
    if(verbose) message("\n\r Compute a distance matrix based on dynamic time wrapping.")
    dat.dist <- dtwDist(mx=dat, method="DTW")
    dat.dist <- as.dist(dat.dist)
  } else if(dist.method == "cor.diss") {
    if(verbose) message("\n\r Compute a distance matrix based on correlation dissimilarity.")
    dat.dist <- TSclust::diss(dat, "COR")
  }

  #### CLUSTER
  if(cluster.method == "hclust") {
    if(verbose) message("\n\r Applying Hierarchical Clustering")
    dat.hc <- hclust(d=dat.dist, method = hclust.algorithm, ...)
    dat.hc.mem <- cutree(dat.hc, k = K)
    centers <- NULL
    for(k in 1:K){
      centers <- rbind(centers, colMeans(dat[dat.hc.mem == k, , drop = FALSE]))
    }
    dat.hc.centers <- hclust(dist(centers), method = hclust.algorithm, members = table(dat.hc.mem))

    if(verbose) {
      opar <- par(mfrow = c(1, 2))
      plot(dat.hc,  labels = FALSE, hang = -1, main = "Hierarchical Tree")
      plot(dat.hc.centers, labels = FALSE, hang = -1, main = "Cluster Centers")
      par(opar)
    }

    output <- list(dat.dist = dat.dist,
                   cluster.obj = dat.hc,
                   membership = dat.hc.mem)
  } else if(cluster.method == "kmeans") {
    if(verbose) message("\n\r Applying Kmeans Clustering")
    if(missing(kmeans.centers)) {
      dat.kmeans <- kmeans(dat, centers = K, algorithm = kmeans.algorithm, ...)
    } else {
      dat.kmeans <- kmeans(dat, algorithm = kmeans.algorithm, centers = kmeans.centers, ...)
    }

    output <- list(cluster.obj= dat.kmeans,
                   membership = dat.kmeans$cluster)
  } else {
    stop("\n\r No Clustering Applied")
  }
  return(output)
}
