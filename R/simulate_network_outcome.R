#' Simulate network-level potential outcomes and observed outcomes under interference
#'
#' This function simulates binary potential outcomes under a network
#' interference setting and then generates observed treatment assignment and
#' observed outcomes from those potential outcomes.
#'
#' For each individual, the potential outcome depends on their own treatment
#' status, the proportion of treated neighbors in the network, and
#' user-specified covariates in the outcome model.
#'
#' The input \code{net_obj} is typically produced by \code{build_network()} and
#' must contain both the network object and the associated node-level data. The
#' node-level degree variable \code{na} is used to define
#' \code{interference_set_size}.
#'
#' If \code{L_i} is included in either \code{outcome_covars} or
#' \code{treatment_covars} but is not already present in \code{net_obj$data},
#' it is generated as a Bernoulli random variable with probability \code{p_L}.
#'
#' The function first constructs a potential outcome matrix for each individual,
#' where rows correspond to individual treatment values (\code{A_i = 0,1}) and
#' columns correspond to the number of treated neighbors
#' (\code{0, 1, \ldots, interference_set_size}).
#'
#' It then simulates treatment assignment from a logistic propensity score
#' model, computes observed network exposure from each node's neighbors, and
#' determines the observed outcome under the realized treatment and neighbor
#' exposure.
#'
#' @param net_obj A list containing the network object and associated node-level
#' data. It must include \code{net0}, an \code{igraph} object, and \code{data},
#' a data frame with one row per node. The data frame must contain \code{na},
#' which gives each node's interference set size (degree). Additional covariates
#' may also be included.
#' @param p_L A numeric scalar giving the Bernoulli probability used to generate
#' \code{L_i} when needed and not already present in \code{net_obj$data}.
#' @param outcome_intercept A numeric scalar for the intercept in the outcome
#' model.
#' @param beta_a A numeric scalar for the main effect of individual treatment in
#' the outcome model.
#' @param beta_abar A numeric scalar for the main effect of the proportion of
#' treated neighbors in the outcome model.
#' @param beta_a_abar A numeric scalar for the interaction between individual
#' treatment and neighbor treatment proportion in the outcome model.
#' @param outcome_covars A character vector giving the names of covariates to be
#' included in the outcome model.
#' @param outcome_coefs A numeric vector of coefficients for
#' \code{outcome_covars}. Must have the same length as \code{outcome_covars}.
#' @param propensity_intercept A numeric scalar for the intercept in the
#' treatment assignment model.
#' @param treatment_covars A character vector giving the names of covariates to
#' be included in the treatment assignment model.
#' @param treatment_coefs A numeric vector of coefficients for
#' \code{treatment_covars}. Must have the same length as
#' \code{treatment_covars}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{PO}}{A list of length \code{nrow(net_obj$data)}. Each element
#'   is a \eqn{2 \times (N_i + 1)} matrix of binary potential outcomes for one
#'   individual. Rows correspond to \code{A_i = 0} and \code{A_i = 1}, and
#'   columns correspond to the number of treated neighbors.}
#'   \item{\code{dat.network}}{A data frame containing the node-level network
#'   data together with:
#'   \describe{
#'     \item{\code{interference_set_size}}{number of neighbors for each node}
#'     \item{\code{L_i}}{generated covariate, if needed}
#'     \item{\code{Ai}}{simulated individual treatment}
#'     \item{\code{treated_neighbors}}{number of treated neighbors in the network}
#'     \item{\code{Abar_Ni}}{proportion of treated neighbors}
#'     \item{\code{Y_obs}}{observed outcome under the realized treatment and neighbor exposure}
#'   }}
#' }
#'
#' The original degree variable \code{na} is removed from the returned network
#' data after being copied to \code{interference_set_size}.
#'
#' @examples
#' net_obj <- build_network()
#' sim_out <- simulate_network_outcome(net_obj)
#'
#' @export
simulate_network_outcome <- function(net_obj,
                                     p_L = 0.4,
                                     outcome_intercept = 0.5,
                                     beta_a = -0.5,
                                     beta_abar = -0.5,
                                     beta_a_abar = 0.25,
                                     outcome_covars = "L_i",
                                     outcome_coefs = 0.25,
                                     propensity_intercept = -0.5,
                                     treatment_covars = "L_i",
                                     treatment_coefs = 1) {

  if (!is.list(net_obj) || is.null(net_obj$net0) || is.null(net_obj$data)) {
    stop("net_obj must be a list containing 'net0' and 'data'.")
  }

  net0 <- net_obj$net0
  dat.network <- net_obj$data

  if (!("na" %in% names(dat.network))) {
    stop("net_obj$data must contain 'na'.")
  }

  dat.network$interference_set_size <- dat.network$na

  if ("L_i" %in% c(outcome_covars, treatment_covars) &&
      !("L_i" %in% names(dat.network))) {
    dat.network$L_i <- rbinom(nrow(dat.network), 1, p_L)
  }

  data <- dat.network
  n <- nrow(data)
  Ni <- data$interference_set_size

  if (length(outcome_covars) != length(outcome_coefs)) {
    stop("length(outcome_covars) must equal length(outcome_coefs)")
  }
  if (length(treatment_covars) != length(treatment_coefs)) {
    stop("length(treatment_covars) must equal length(treatment_coefs)")
  }

  if (!all(outcome_covars %in% names(data))) {
    stop("Some outcome_covars are not in dat.network")
  }
  if (!all(treatment_covars %in% names(data))) {
    stop("Some treatment_covars are not in dat.network")
  }

  PO <- vector("list", length = n)
  for (i in seq_len(n)) {
    PO[[i]] <- matrix(NA_integer_, nrow = 2, ncol = Ni[i] + 1)
  }

  for (i in seq_len(n)) {
    X_outcome_i <- as.numeric(data[i, outcome_covars, drop = TRUE])

    for (j in 0:1) {
      for (k in 0:Ni[i]) {
        abar <- if (Ni[i] == 0) 0 else k / Ni[i]

        eta <- outcome_intercept +
          beta_a * j +
          beta_abar * abar +
          beta_a_abar * j * abar +
          sum(outcome_coefs * X_outcome_i)

        p <- plogis(eta)
        PO[[i]][j + 1, k + 1] <- rbinom(1, 1, p)
      }
    }
  }

  X_treat <- as.matrix(data[, treatment_covars, drop = FALSE])
  pA <- plogis(propensity_intercept + X_treat %*% treatment_coefs)
  dat.network$Ai <- rbinom(nrow(dat.network), size = 1, prob = as.numeric(pA))

  dat.network$treated_neighbors <- vapply(
    seq_len(nrow(dat.network)),
    function(i) {
      nbrs <- igraph::neighbors(net0, i)
      sum(dat.network$Ai[as.integer(nbrs)])
    },
    numeric(1)
  )

  dat.network$Abar_Ni <- ifelse(
    dat.network$interference_set_size == 0,
    0,
    dat.network$treated_neighbors / dat.network$interference_set_size
  )

  dat.network$Y_obs <- vapply(
    seq_len(nrow(dat.network)),
    function(i) {
      PO_i <- PO[[i]]
      PO_i[
        dat.network$Ai[i] + 1,
        dat.network$treated_neighbors[i] + 1
      ]
    },
    numeric(1)
  )

  dat.network$na <- NULL

  return(list(
    PO = PO,
    dat.network = dat.network
  ))
}
