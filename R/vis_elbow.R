#' Visualize the percent variance explained by PCs
#'
#' Plot the classic elbow plot, showing the percent variance explained by PCs.
#'
#' The input data is either centered and/or scaled.
#' SVD/PCA is computed and the percent variacnes explained by the top \code{r} PCs are plotted.
#' If the denoising step is used, please input the denoised data.
#'
#' @param dat a data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param center.dat a logical specifying to center the data. By default, \code{TRUE}.
#' @param scale.dat a logical specifying to scale the data. By default, \code{FALSE}.
#' @param r the top \code{r} PCs to display. by default, showing all PCs.
#' @param pt.size the size of data points.
#' @param title the title for the resulting plot.
#'
#' @return \code{pca.elbow} returns a ggplot2 object.
#'
#' @importFrom corpcor fast.svd
#' @export pca.elbow
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @examples
#' \dontrun{
#' data(cys_optm)
#' meta <- cys_optm[,1:4]
#' optm <- log(cys_optm[meta$Select,5:10])
#' optm <- t(scale(t(optm), scale=TRUE, center=TRUE))
#' days <- as.numeric(colnames(optm))
#'
#' pca.elbow(dat=optm,
#' center.dat=TRUE,
#' scale.dat=FALSE,
#' title="Percent variance explained by 10 PCs")
#'}
pca.elbow <- function (dat,
                       center.dat=TRUE,
                       scale.dat=FALSE,
                       r=NULL,
                       pt.size = 1,
                       title="")
{
  m <- nrow(dat)
  n <- ncol(dat)
  if(scale.dat | center.dat) dat <- t(scale(t(dat),center=center.dat,scale=scale.dat))
  svd.obj <- fast.svd(dat)
  if(is.null(r)) r = length(svd.obj$d)

  pve <- svd.obj$d[1:r]^2 / sum(svd.obj$d^2) * 100
  data.plot <- data.frame(PVE=pve, PC=1:length(svd.obj$d[1:r]), pt.size=pt.size)
  p <- ggplot(data = data.plot, mapping = aes(x = PC, y = PVE)) +
    geom_point(size = pt.size) + theme_bw() +
    theme(legend.title = element_blank()) + ggtitle(title)
  return(p)
}

#' Visually inspect the optimal numbers of clusters
#'
#' A convenient wrapper for \code{factoextra::fviz_nbclust}, which plot within cluster sums of squares, average silhouette and gap statistics for clustering.
#' The required parameters and their corresponding documentations are simplified versions of that in  \code{factoextra::fviz_nbclust}.
#'
#' @param dat an input data.
#' @param FUNcluster a partitioning function which accepts as first argument a (data) matrix like x, second argument, say k, k >= 2, the number of clusters desired, and returns a list with a component named cluster which contains the grouping of observations. By default, "kmeans".
#' @param method the method to be used for estimating the optimal number of clusters. Possible values are "silhouette" (for average silhouette width), "wss" (for total within sum of square) and "gap_stat" (for gap statistics).
#' @param k.max the maximum number of clusters to consider. By default, 10.
#' @param ... optional arguments. See or directly use \code{factoextra::fviz_nbclust}.
#'
#' @return \code{cluster.elbow} returns a ggplot2 object.
#'
#' @importFrom factoextra fviz_nbclust
#' @export cluster.elbow
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @examples
#' \dontrun{
#' data(coptm)
#' cluster.elbow(dat=optm, FUNcluster=kmeans, method="wss", k.max=10, linecolor="black")
#'}
cluster.elbow <- function (dat,
                           FUNcluster,
                           method = c("wss", "silhouette", "gap_stat"),
                           k.max = 10,
                           ...)
{
  p <- factoextra::fviz_nbclust(x=dat, FUNcluster=FUNcluster, method=method, k.max=k.max, ...)
  return(p)
}
