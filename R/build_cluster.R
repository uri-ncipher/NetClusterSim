#' Generate Clustered Data with Interference Set Sizes
#'
#' This function generates clustered data where each subject belongs to a cluster,
#' and each cluster has a size drawn from a specified distribution (default: Poisson).
#' The output includes subject IDs, cluster IDs, and interference set sizes.
#'
#' Clusters smaller than \code{min_size} are removed.
#'
#' @param K Integer. Number of initial clusters to generate.
#' @param dist Character. Distribution for cluster sizes.
#'   Options: \code{"poisson"} (default), \code{"nbinom"}, \code{"uniform"}, \code{"fixed"}.
#' @param lambda Numeric. Mean parameter for Poisson distribution (default = 10).
#' @param size Numeric. Size parameter for negative binomial distribution.
#' @param prob Numeric. Success probability for negative binomial distribution.
#' @param min_size Integer. Minimum cluster size to retain (default = 2).
#' @param max_size Integer. Maximum cluster size (used for uniform distribution).
#' @param fixed_size Integer. Fixed cluster size (used when \code{dist = "fixed"}).
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{subject_id}{Unique subject identifier}
#'   \item{cluster_id}{Cluster membership ID}
#'   \item{interference_set_size}{Cluster size minus one}
#' }
#'
#' Attributes:
#' \describe{
#'   \item{cluster_sizes}{Vector of retained cluster sizes}
#'   \item{dist}{Distribution used}
#'   \item{n_clusters_kept}{Number of clusters after filtering}
#'   \item{original_K}{Original number of clusters}
#' }
#'
#' @examples
#' # Default (Poisson clusters)
#' dat <- build_cluster(K = 1000, lambda = 10, seed = 123)
#' head(dat)
#'
#' # Negative binomial clusters
#' dat_nb <- build_cluster(K = 1000, dist = "nbinom",
#'                                 size = 5, prob = 0.4, seed = 123)
#'
#' # Uniform cluster sizes
#' dat_unif <- build_cluster(K = 1000, dist = "uniform",
#'                                   min_size = 2, max_size = 15, seed = 123)
#'
#' # Fixed cluster size
#' dat_fixed <- build_cluster(K = 1000, dist = "fixed",
#'                                    fixed_size = 8, seed = 123)
#'
#' @export



build_cluster <- function(K = 1000,
                                  dist = c("poisson", "nbinom", "uniform", "fixed"),
                                  lambda = 10,
                                  size = NULL,
                                  prob = NULL,
                                  min_size = 2,
                                  max_size = NULL,
                                  fixed_size = NULL,
                                  seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  dist <- match.arg(dist)

  if (!is.numeric(K) || length(K) != 1 || K <= 0) {
    stop("K must be a single positive number.")
  }
  if (!is.numeric(min_size) || length(min_size) != 1 || min_size < 1) {
    stop("min_size must be a single number >= 1.")
  }

  K <- as.integer(K)
  min_size <- as.integer(min_size)

  cluster_size <- switch(
    dist,
    poisson = {
      if (is.null(lambda) || lambda <= 0) stop("lambda must be > 0 for poisson.")
      stats::rpois(K, lambda)
    },
    nbinom = {
      if (is.null(size) || size <= 0) stop("size must be > 0 for nbinom.")
      if (is.null(prob) || prob <= 0 || prob >= 1) stop("prob must be in (0, 1) for nbinom.")
      stats::rnbinom(K, size = size, prob = prob)
    },
    uniform = {
      if (is.null(max_size) || max_size < min_size) {
        stop("max_size must be >= min_size for uniform.")
      }
      sample(min_size:max_size, K, replace = TRUE)
    },
    fixed = {
      if (is.null(fixed_size) || fixed_size < 1) {
        stop("fixed_size must be >= 1 for fixed.")
      }
      rep(as.integer(fixed_size), K)
    }
  )

  keep <- which(cluster_size >= min_size)

  if (length(keep) == 0) {
    stop("No clusters satisfy the minimum size condition.")
  }

  kept_sizes <- cluster_size[keep]

  dat.cluster <- data.frame(
    cluster_id = rep(seq_along(keep), times = kept_sizes)
  )

  dat.cluster$subject_id <- seq_len(nrow(dat.cluster))
  dat.cluster$interference_set_size <- kept_sizes[dat.cluster$cluster_id] - 1

  dat.cluster <- dat.cluster[, c("subject_id", "cluster_id", "interference_set_size")]

  attr(dat.cluster, "cluster_sizes") <- kept_sizes
  attr(dat.cluster, "dist") <- dist
  attr(dat.cluster, "n_clusters_kept") <- length(keep)
  attr(dat.cluster, "original_K") <- K

  return(dat.cluster)
}

