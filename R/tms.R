#' Identify the temporal molecular signatures
#'
#' Denoise, classify, and evaluate variables (biomarkers) from time course data such as proteomics and other high-throughput technologies.
#'
#' This function combines multiple steps. For more options and fine-tuning, please use individual functions in `tms` package.
#' \code{tms} attempts to identify temporal dynamics, by clustering denoise and/or time-wrapped data.
#' This requires the user to input the data (\code{dat}) where the rows and columns are variables (e.g., genes, proteins) and observations taken at different time points, respectively.
#' Correspondingly, \code{timepoints} is a vector of actual time points (e.g., hours, days) corresponding to the columns of \code{dat}.
#'
#' This function goes through the following steps:
#' \itemize{
#' \item denoise temporal data using cubic splines (\code{denoise_spline}) or PCA (\code{denoise_pca}).
#' When \code{denoise="cubic.spline"} is chosen, individual degrees of freedom can be chosen by cross validation by setting \code{denoise.parameter = "cv"}. If there should only one tuning parameter for all variables, set \code{denoise.parameter = "cv.global"}.
#' \item cluster variables (proteins, genes, etc) using \code{cluster.method} based on \code{dist.method} and \code{K} clusters.
#' When using dynamic time wrapping (DTW) \code{dist.method = "dtw"}, hierachical clustering is applied.
#' K-means clustering (\code{dist.method = "kmeans"}) does not return a distance matrix.
#' \item evaluate the cluster memberships of variables (e.g., proteins or genes) by the jackstraw tests.
#' The jackstraw returns p-values and posterior probabilities that variables should be included in their given clusters.
#' }
#'
#' This work is motivated by identifying reliable molecular signatures from time-series proteomics data of optm occupancies in the cardiovascular mouse model (see Wang et al. (2018))
#' Last but not least, modeling and classifying high-dimensional temporal data is notoriously challenging. This package aim to provide an analysis pipeline that is relatively robust and non-parametric, while accounting for typical -omic study involving complex phenotypes. For further implementations of related methods, see \code{TSclust} and \code{TSdist}.
#'
#' @param dat a data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param timepoints a vector of time points for columns of dat.
#' @param center.dat a logical specifying to center the input and denoised data. By default, \code{TRUE}.
#' @param scale.dat a logical specifying to scale the input and denoised data. By default, \code{FALSE}.
#' @param denoise a denoising method. By default, fitting a cubic spline.
#' @param denoise.parameter a parameter for a denoising method, such as the degree of freedom in spline.smooth, the number of significant PCs in PCA.
#' @param dist.method a distance method for time course data, resulting in a \code{m * m} distance matrix for rows. 'dtw' for dynamic time wrapping or 'cor.diss' for correlation-based dissimilarities.
#' @param cluster.method a clustering method.
#' @param K a number of clusters.
#' @param evaluate a logical specifying to evaluate the cluster membership with the jackstraw tests. By default, \code{FALSE}.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param ... optional arguments.
#'
#' @return \code{tms} returns a list consisting of
#' \item{denoised}{\code{m * n} denoised data.}
#' \item{dat.dist}{\code{m * m} distance matrix. Only returns when using hclust}
#' \item{cluster.obj}{an object returned from clustering the denoised data.}
#' \item{membership}{a vector of length \code{m}, identities of clusters.}
#' \item{evaluated}{an object returned from applying the jackstraw tests for clusters.}
#'
#' @references Identifying temporal molecular signatures underlying cardiovascular diseases. In preparation.
#' @references J Wang, H Choi, NC Chung, Q Cao, DCM Ng, B Mirza, SB Scruggs, D Wang, AO Garlid, P Ping (2018). Integrated dissection of the cysteine oxidative post-translational modification proteome during cardiac hypertrophy. Journal of Proteome Research.
#' @references NC Chung (2020). Statistical significance of cluster membership for unsupervised evaluation of single cell identities. Bioinformatics
#'
#' @export tms
#' @import splines
#' @importFrom jackstraw jackstraw_kmeans
#' @importFrom corpcor fast.svd
#' @importFrom dtw dtwDist
#' @importFrom TSclust diss
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'@examples
#'\dontrun{
#' data(cys_optm)
#' meta <- cys_optm[,1:4]
#' optm <- log(cys_optm[meta$Select,5:10])
#' optm <- t(scale(t(optm), scale=TRUE, center=TRUE))
#' days <- as.numeric(colnames(optm))
#'
#' optm.tms <- tms(optm,
#'                 timepoints = days,
#'                 center.dat = TRUE,
#'                 scale.dat = TRUE,
#'                 denoise = c("smooth.spline"),
#'                 denoise.parameter=c("cv"),
#'                 dist.method = "cor.diss",
#'                 cluster.method = c("kmeans"),
#'                 K = 5,
#'                 evaluate = TRUE,
#'                 verbose = TRUE,
#'                 seed = 1
#'                )
#'
#' # see the elbow plot
#' cluster.elbow(dat=optm.tms$denoised, FUNcluster=kmeans, method="wss", k.max=10, linecolor="black")
#'
#' # make the cluster figure
#' optm.fig <- vis_cluster(optm.tms$denoised, group=optm.tms$membership)
#'
#' # to modify/polish the figure (ggplot2 object)
#' optm.fig <- optm.fig + labs(y="Log-transformed Occupancy Ratio", x="Time (day)", title="All O-PTMs") + ylim(-2,2) + facet_wrap(~ cluster,nrow=1,ncol=6)
#'
#' # filter the data based on jackstraw PIP and make a figure
#' library(jackstraw)
#' optm.pip <- pip(optm.tms$evaluated$p.F, pi0=sum(optm.tms$evaluated$p.F > .05)/length(optm.tms$evaluated$p.F))
#' hist(optm.pip,100,col="black")
#' optm.pip.fig <- vis_cluster(optm.tms$denoised[optm.pip > .9,], group=optm.tms$membership[optm.pip > .9])
#' optm.pip.fig <- optm.pip.fig + labs(y="Log-transformed Occupancy Ratio", x="Time (day)", title="O-PTMs with PIP > 0.9") + ylim(-2,2) + facet_wrap(~ cluster,nrow=1,ncol=6)
#'
#' library(cowplot)
#' optm.fig / optm.pip.fig
#'}
tms <- function(dat,
                timepoints=NULL,
                center.dat = TRUE,
                scale.dat = FALSE,
                denoise = c("smooth.spline", "pca", "none"),
                denoise.parameter = c("cv", "cv.global"),
                dist.method = c("euclidean", "cor.diss", "dtw"),
                cluster.method = c("kmeans", "hclust"),
                K,
                evaluate = TRUE,
                verbose = FALSE,
                seed = NULL) {

  if (is.null(seed)) set.seed(seed)
  if(verbose)  {
    message("\n\r The rows are variables (e.g., proteins) and the columns are time points.")
  }
  m <- nrow(dat)
  n <- ncol(dat)
  if(scale.dat | center.dat) dat <- t(scale(t(dat),center=center.dat,scale=scale.dat))

  if(cluster.method == "kmeans") {
    if(verbose) message("You have chosen cluster.method == 'kmeans'. Therefore, it won't calculate a distance matrix.")
    dist.method == "none"
  }
  # sanity check
  if(n != length(timepoints)) stop("The number of time points must match the number of columns in the input data.")

  #### 1. DENOISE
  if(denoise == "smooth.spline") {
    dat.denoised <- denoise_spline(dat=dat, timepoints=timepoints, dof=denoise.parameter,
                                  center.dat = center.dat, scale.dat = scale.dat,
                                  verbose = verbose, seed = NULL)
  } else if(denoise == "pca"){
    dat.denoised <- denoise_pca(dat=dat, timepoints=timepoints, r=denoise.parameter,
                               center.dat = center.dat, scale.dat = scale.dat,
                               verbose = verbose, seed = NULL)
  } else {
    stop("No denoising selected. If you would like to proceed without denoising, simply use individual functions such as `tms::cluster`.")
  }

  if(scale.dat | center.dat) dat.denoised <- t(scale(t(dat.denoised),center=center.dat,scale=scale.dat))

  #### 2. CLUSTER
  dat.clustered <- cluster(dat.denoised,
                           timepoints=timepoints,
                           K=K,
                           dist.method = dist.method,
                           cluster.method = cluster.method,
                           center.dat = center.dat,
                           scale.dat = scale.dat,
                           verbose = verbose)

  #### 3. EVALUATE
  if(evaluate) {
    if(cluster.method == "kmeans") {
      message("Evaluate the cluster memberships with the jackstraw tests")
      dat.evaluated <- jackstraw_kmeans(dat = dat.denoised,
                                        kmeans.dat = dat.clustered$cluster.obj,
                                        center = center.dat,
                                        verbose = FALSE)
      dat.evaluated[1] <- NULL #drop the match.call() from the jackstraw for a cleaner output
    } else {
      stop("Currently the evaluation step is only available for kmeans clustering.")
    }
  }

  #### RETURN
  output <- append(dat.clustered,
                   list(denoised=dat.denoised,
                        evaluated=dat.evaluated)
                   )
  return(output)
}
