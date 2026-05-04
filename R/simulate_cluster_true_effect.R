#' Compute true causal effects under clustered interference
#'
#' This function computes true causal estimands under a clustered interference
#' setting using the output from \code{simulate_cluster_outcome()}.
#'
#' The input must contain a list of individual potential outcome matrices
#' together with the corresponding clustered data structure. For each
#' individual, the potential outcome matrix is assumed to have two rows,
#' corresponding to individual treatment values \code{A_i = 0} and
#' \code{A_i = 1}, and \code{N_i + 1} columns, corresponding to the number of
#' treated neighbors in the individual's interference set. The function averages
#' over the distribution of treated neighbors under user-specified treatment
#' allocation probabilities.
#'
#' The following causal estimands are computed:
#' \describe{
#'   \item{\code{DE}}{Direct effect at a fixed allocation \eqn{\alpha}}
#'   \item{\code{IE}}{Indirect effect comparing two allocations}
#'   \item{\code{TE}}{Total effect comparing treatment and allocation jointly}
#'   \item{\code{OE}}{Overall effect comparing two marginal allocations}
#' }
#'
#' Cluster-level means are computed by first averaging within each cluster and
#' then averaging across clusters.
#'
#' @param sim_out A list returned by \code{simulate_cluster_outcome()}.
#' It must contain:
#' \describe{
#'   \item{\code{PO}}{a list of potential outcome matrices, one for each individual}
#'   \item{\code{dat.cluster}}{a data frame containing the clustered structure,
#'   including \code{subject_id}, \code{cluster_id}, and
#'   \code{interference_set_size}}
#' }
#' @param allocations A numeric vector of treatment allocation probabilities.
#' Each value must lie in \eqn{[0,1]}. Defaults to
#' \code{c(0.25, 0.5, 0.75)}.
#'
#' @return A list with class \code{"cluster_true_effect"} containing:
#' \describe{
#'   \item{\code{allocations}}{the input allocation probabilities}
#'   \item{\code{DE}}{a data frame of true direct effects for each allocation}
#'   \item{\code{IE}}{a data frame of true indirect effects for each pair of allocations}
#'   \item{\code{TE}}{a data frame of true total effects for each pair of allocations}
#'   \item{\code{OE}}{a data frame of true overall effects for each pair of allocations}
#' }
#'
#' The \code{DE} data frame contains:
#' \describe{
#'   \item{\code{effect}}{effect label, always \code{"DE"}}
#'   \item{\code{alpha}}{allocation probability}
#'   \item{\code{true_effect}}{true direct effect at that allocation}
#' }
#'
#' The \code{IE}, \code{TE}, and \code{OE} data frames contain:
#' \describe{
#'   \item{\code{effect}}{effect label}
#'   \item{\code{alpha1}}{first allocation probability}
#'   \item{\code{alpha0}}{second allocation probability}
#'   \item{\code{true_effect}}{true causal effect comparing the two allocations}
#' }
#'
#' @examples
#' dat.cluster <- build_cluster()
#' sim_out <- simulate_cluster_outcome(dat.cluster)
#' simulate_cluster_true_effect(sim_out)
#'
#' @export

simulate_cluster_true_effect <- function(sim_out,
                                         allocations = c(0.25, 0.5, 0.75)) {

  if (!is.list(sim_out) || is.null(sim_out$PO) || is.null(sim_out$dat.cluster)) {
    stop("sim_out must be a list containing 'PO' and 'dat.cluster'.")
  }

  PO <- sim_out$PO
  dat.cluster <- sim_out$dat.cluster

  req_cols <- c("subject_id", "cluster_id", "interference_set_size")
  if (!all(req_cols %in% names(dat.cluster))) {
    stop("dat.cluster must contain 'subject_id', 'cluster_id', and 'interference_set_size'.")
  }

  if (!is.list(PO)) {
    stop("PO must be a list of potential outcome matrices.")
  }

  if (length(PO) != nrow(dat.cluster)) {
    stop("length(PO) must equal nrow(dat.cluster).")
  }

  if (length(allocations) < 1) {
    stop("allocations must contain at least one value.")
  }

  if (any(allocations < 0 | allocations > 1)) {
    stop("allocations must be between 0 and 1.")
  }

  y_bar_ind <- function(PO, i, a, alpha) {
    M <- PO[[i]]
    n_nei <- ncol(M) - 1
    Arow <- a + 1

    out <- 0
    for (k in 0:n_nei) {
      out <- out + M[Arow, k + 1] * stats::dbinom(k, size = n_nei, prob = alpha)
    }
    out
  }

  y_bar <- function(PO, dat.cluster, a, alpha) {
    cl_ids <- unique(dat.cluster$cluster_id)
    out <- 0

    for (cl in cl_ids) {
      idx <- which(dat.cluster$cluster_id == cl)
      cl_mean <- mean(vapply(idx, function(i) y_bar_ind(PO, i, a, alpha), numeric(1)))
      out <- out + cl_mean
    }

    out / length(cl_ids)
  }

  y_bar_ind_margin <- function(PO, i, alpha) {
    M <- PO[[i]]
    n_nei <- ncol(M) - 1

    out <- 0
    for (a in 0:1) {
      pa <- if (a == 1) alpha else (1 - alpha)
      for (k in 0:n_nei) {
        out <- out + M[a + 1, k + 1] *
          pa *
          stats::dbinom(k, size = n_nei, prob = alpha)
      }
    }
    out
  }

  y_bar_margin <- function(PO, dat.cluster, alpha) {
    cl_ids <- unique(dat.cluster$cluster_id)
    out <- 0

    for (cl in cl_ids) {
      idx <- which(dat.cluster$cluster_id == cl)
      cl_mean <- mean(vapply(idx, function(i) y_bar_ind_margin(PO, i, alpha), numeric(1)))
      out <- out + cl_mean
    }

    out / length(cl_ids)
  }

  DE_true <- function(alpha) {
    y_bar(PO, dat.cluster, a = 1, alpha = alpha) -
      y_bar(PO, dat.cluster, a = 0, alpha = alpha)
  }

  IE_true <- function(alpha1, alpha0) {
    y_bar(PO, dat.cluster, a = 0, alpha = alpha1) -
      y_bar(PO, dat.cluster, a = 0, alpha = alpha0)
  }

  TE_true <- function(alpha1, alpha0) {
    y_bar(PO, dat.cluster, a = 1, alpha = alpha1) -
      y_bar(PO, dat.cluster, a = 0, alpha = alpha0)
  }

  OE_true <- function(alpha1, alpha0) {
    y_bar_margin(PO, dat.cluster, alpha = alpha1) -
      y_bar_margin(PO, dat.cluster, alpha = alpha0)
  }

  # true direct effects at each allocation
  de_df <- data.frame(
    effect = "DE",
    alpha = allocations,
    true_effect = vapply(allocations, DE_true, numeric(1))
  )

  # pairwise true effects between allocations
  pair_mat <- utils::combn(allocations, 2)

  ie_df <- data.frame(
    effect = "IE",
    alpha1 = pair_mat[1, ],
    alpha0 = pair_mat[2, ],
    true_effect = vapply(seq_len(ncol(pair_mat)), function(j) {
      IE_true(pair_mat[1, j], pair_mat[2, j])
    }, numeric(1))
  )

  te_df <- data.frame(
    effect = "TE",
    alpha1 = pair_mat[1, ],
    alpha0 = pair_mat[2, ],
    true_effect = vapply(seq_len(ncol(pair_mat)), function(j) {
      TE_true(pair_mat[1, j], pair_mat[2, j])
    }, numeric(1))
  )

  oe_df <- data.frame(
    effect = "OE",
    alpha1 = pair_mat[1, ],
    alpha0 = pair_mat[2, ],
    true_effect = vapply(seq_len(ncol(pair_mat)), function(j) {
      OE_true(pair_mat[1, j], pair_mat[2, j])
    }, numeric(1))
  )

  out <- list(
    allocations = allocations,
    DE = de_df,
    IE = ie_df,
    TE = te_df,
    OE = oe_df
  )

  class(out) <- "cluster_true_effect"
  out
}
