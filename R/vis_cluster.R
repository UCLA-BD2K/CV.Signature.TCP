#' Visualize the grouped projectories
#'
#' Plot the time course projectories that are grouped by clustering or user-defined groups.
#'
#' @param dat a time-series data matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param group a vector defining the groups.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param seed a seed for the random number generator.
#' @param center.dat a logical specifying to center the input and denoised data. By default, \code{TRUE}.
#' @param scale.dat a logical specifying to scale the input and denoised data. By default, \code{FALSE}.
#' @param ... optional arguments.
#'
#' @return \code{vis_cluster} returns a ggplot2 object.
#'
#' @import ggplot2
#' @import magrittr
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom tibble add_column
#' @export vis_cluster
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
vis_cluster <- function(dat,
                      group = NULL,
                      verbose = TRUE,
                      seed = NULL,
                      center.dat = TRUE,
                      scale.dat = FALSE,
                      ...) {

  if (is.null(seed)) set.seed(seed)
  if (is.null(group)) stop("Groups of rows must be supplied as a vector")
  m <- nrow(dat)
  n <- ncol(dat)

  if(scale.dat | center.dat) dat <- t(scale(t(dat),center=center.dat,scale=scale.dat))
  if(length(group) != nrow(dat)) stop("The length of group must equal the number of rows in data.")

  K <- length(unique(group))
  # PLOT ALL DATA POINTS WITH CLUSTER CENTERS
  dat.plot <- dat
  id <- dat.plot %>% as.data.frame() %>% tibble::rownames_to_column() %>%
    add_column(cluster = group)
  id <- melt(id, variable.name="day", value.name="ratio", id.vars=c("rowname","cluster"))
  id$day <- as.numeric(as.character(id$day))
  id$rowname <- as.factor(id$rowname)

  gd <- id %>%
    group_by(cluster, day) %>%
    summarise(ratio_mean = mean(ratio))

  g <- ggplot(id, aes(x = day, y = ratio)) +
    geom_point(alpha = .3, size=.5) +
    geom_line(aes(group=rowname), alpha = .3) +
    geom_line(data = gd, aes(x = day, y = ratio_mean, group=1), alpha = 1, size = 1, col="darkorange", linetype="solid") +
    theme_bw() +
    facet_wrap(~ cluster)

  return(g)
}

#' Visualize the projectories of two datasts
#'
#' Plot two sets of projectories next to each other, grouped by clusters.
#'
#' @param x a first matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param x.cluster a clustering object from the first data matrix.
#' @param y a second matrix with \code{m} biomarkers as rows, over \code{n} time points (columns).
#' @param y.cluster a clustering object from the second data matrix.
#' @param verbose a logical specifying to print the computational progress. By default, \code{FALSE}.
#' @param print.plot a logical specifying to print the resulting plot onto the screen. By default, \code{FALSE}.
#' @param return.separate a logical specifying to return two ggplot2 objects separately. By default, \code{FALSE}.
#' @param ... optional arguments.
#'
#' @return \code{vis_cluster_multi} returns a ggplot2 object.
#'
#' @importFrom dplyr group_by
#' @importFrom tibble add_column
#' @import patchwork
#' @import ggplot2
#' @export vis_cluster_multi
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
vis_cluster_multi <- function(x, x.cluster,
                         y, y.cluster,
                         verbose = TRUE,
                         print.plot = FALSE,
                         return.separate = FALSE,
                         ...) {

  m <- nrow(x)
  n <- ncol(y)

  # ensure the inputs are valid
  if(nrow(x.cluster$centers) != nrow(y.cluster$centers)) {
    stop("The number of cluster centers must be identical.")
  }
  if(length(unique(x.cluster$cluster)) != length(unique(y.cluster$cluster))) {
    message("The number of unique cluster assignments are not identical. Make sure you know what you are doing.")
  }
  if(!all(dim(x) == dim(y))) {
    stop("The dimensions of the data x and the data y must be identical.")
  }
  if(!all(rownames(x) == rownames(y))) {
    stop("The rownames of the data x and the data y must be identical.")
  }
  if(!all(colnames(x) == colnames(y))) {
    stop("The colnames of the data x and the data y must be identical.")
  }

  # PLOT x
  xi <- x %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    add_column(cluster = x.cluster$cluster)
  xi <- melt(xi, variable.name="time", value.name="value", id.vars=c("rowname","cluster"))
  xi$time <- as.numeric(as.character(xi$time))

  xg <- xi %>%
    group_by(cluster, time) %>%
    summarise(value_mean = mean(value))

  plotx <- ggplot(xi, aes(x = time, y = value)) +
    geom_point(alpha = .3, size=.5) +
    geom_line(aes(group=rowname), alpha = .3) +
    geom_line(data = xg, aes(x = time, y = value_mean, group=1), alpha = 1, size = 1, col="darkorange", linetype="solid") +
    theme_bw() +
    facet_wrap(~ cluster, scales ="free_y", ncol = 1) +
    ggtitle("x")

  # PLOT y
  yi <- y %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    add_column(cluster = y.cluster$cluster)
  yi <- melt(yi, variable.name="time", value.name="value", id.vars=c("rowname","cluster"))
  yi$time <- as.numeric(as.character(yi$time))

  yg <- yi %>%
    group_by(cluster, time) %>%
    summarise(value_mean = mean(value))

  ploty <- ggplot(yi, aes(x = time, y = value)) +
    geom_point(alpha = .3, size=.5) +
    geom_line(aes(group=rowname), alpha = .3) +
    geom_line(data = yg, aes(x = time, y = value_mean, group=1), alpha = 1, size = 1, col="darkorange", linetype="solid") +
    theme_bw() +
    facet_wrap(~ cluster, scales="free_y", ncol = 1) +
    ggtitle("y")

  if(print.plot) plotx + ploty
  if(return.separate) {
    return(list(plotx = plotx, ploty = ploty))
  } else {
    return(plotx + ploty)
  }
}
