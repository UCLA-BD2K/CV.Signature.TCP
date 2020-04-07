#' Visualize two PCs
#'
#' Plot two PCs where each rows (biomarkers) are colored by their group
#'
#' @param svd.dat a SVD object. See corpcor::fastsvd.
#' @param ident identities of rows (e.g., clusters).
#' @param dim.1 which dimension of PCA to display as the x-axis
#' @param dim.2 which dimension of PCA to display as the y-axis
#' @param pt.size the size of data points
#' @param pip posterior inclusion probabilities
#' @param no.axes show axes. by default, FALSE
#' @param coord_fixed fix the coordinate. by default, FALSE
#' @param title provide a title
#'
#' @return \code{svd.vis} returns a ggplot2 object.
#'
#' @importFrom dplyr group_by
#' @importFrom tibble add_column
#' @export svd.vis
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
svd.vis <- function (svd.dat, ident, dim.1 = 1, dim.2 = 2, pt.size = 1, pip = 1, no.axes=FALSE, coord_fixed=FALSE, title="")
{
  dim.codes <- paste0("PC", c(dim.1, dim.2))
  data.plot <- data.frame(dim.1=svd.dat$u[,dim.1],dim.2=svd.dat$u[,dim.2])
  data.plot$ident <- ident
  data.plot$pip <- pip
  data.plot$pt.size <- pt.size
  p <- ggplot(data = data.plot, mapping = aes(x = dim.1, y = dim.2))
  if(all(pip==1)) p <- p + geom_point(mapping = aes(colour = factor(x = ident)), size = pt.size)
  if(!all(pip==1)) p <- p + geom_point(mapping = aes(colour = factor(x = ident), alpha = pip), size = pt.size)
  if(coord_fixed) p <- p + coord_fixed()
  p <- p + xlab(label = dim.codes[1]) + ylab(label = dim.codes[2])
  p <- p + theme_bw()
  p <- p + theme(legend.title = element_blank()) + ggtitle(title) + theme(legend.position="bottom")
  if(no.axes) {
    p <- p + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank(),
                   axis.title.x = element_blank(), axis.title.y = element_blank(),
                   panel.background = element_blank(), panel.border = element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   plot.background = element_blank())
  }
  return(p)
}

#' Visualize two dimensions of t-SNE
#'
#' Plot two dimensions of t-SNE where each rows (biomarkers) are colored by their identities
#'
#' @param tsne.obj a t-SNE object.
#' @param ident identities of rows (e.g., clusters).
#' @param dim.1 which dimension of t-SNE to display as the x-axis
#' @param dim.2 which dimension of t-SNE to display as the y-axis
#' @param pt.size the size of data points
#' @param pip posterior inclusion probabilities
#' @param no.axes show axes. by default, FALSE
#' @param coord_fixed fix the coordinate. by default, FALSE
#' @param title provide a title
#'
#' @return \code{tsne.vis} returns a ggplot2 object.
#'
#' @importFrom dplyr group_by
#' @importFrom tibble add_column
#' @export tsne.vis
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
tsne.vis <- function (tsne.obj, ident, dim.1 = 1, dim.2 = 2, pt.size = 1, pip = 1, no.axes=FALSE, coord_fixed=FALSE, title="")
{
  dim.codes <- paste0("t-SNE ", c(dim.1, dim.2))
  data.plot <- data.frame(dim.1=tsne.obj[,dim.1],dim.2=tsne.obj[,dim.2])
  data.plot$ident <- ident
  data.plot$pip <- pip
  data.plot$pt.size <- pt.size
  p <- ggplot(data = data.plot, mapping = aes(x = dim.1, y = dim.2))
  if(all(pip==1)) p <- p + geom_point(mapping = aes(colour = factor(x = ident)), size = pt.size)
  if(!all(pip==1)) p <- p + geom_point(mapping = aes(colour = factor(x = ident), alpha = pip), size = pt.size)
  if(coord_fixed) p <- p + coord_fixed()
  p <- p + xlab(label = dim.codes[1]) + ylab(label = dim.codes[2])
  p <- p + theme_bw() +
    theme(legend.title = element_blank()) + ggtitle(title) +
    theme(legend.position="bottom")
  if(no.axes) {
    p <- p + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank(),
                   axis.title.x = element_blank(), axis.title.y = element_blank(),
                   panel.background = element_blank(), panel.border = element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   plot.background = element_blank())
  }
  return(p)
}
