#' Simulate cluster-level potential outcomes and observed outcomes under interference
#'
#' This function simulates binary potential outcomes under a clustered
#' interference setting and then generates observed treatment assignment and
#' observed outcomes from those potential outcomes.
#'
#' For each individual, the potential outcome depends on their own treatment
#' status, the proportion of treated neighbors within the same cluster, and
#' user-specified covariates in the outcome model.
#'
#' If \code{L_i} is included in either \code{outcome_covars} or
#' \code{treatment_covars} but is not already present in \code{dat.cluster},
#' it is generated as a Bernoulli random variable with probability \code{p_L}.
#'
#' The function first constructs a potential outcome matrix for each individual,
#' where rows correspond to individual treatment values (\code{A_i = 0,1}) and
#' columns correspond to the number of treated neighbors
#' (\code{0, 1, \ldots, interference\_set\_size}).
#'
#' It then simulates treatment assignment from a logistic propensity score
#' model, computes observed cluster exposure, and determines the observed
#' outcome under the realized treatment and neighbor exposure.
#'
#' @param dat.cluster A data frame containing clustered data. It must include
#' at least \code{cluster_id} and \code{interference_set_size}. Additional
#' covariates may also be included.
#' @param p_L A numeric scalar giving the Bernoulli probability used to generate
#' \code{L_i} when needed and not already present in \code{dat.cluster}.
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
#'   \item{\code{PO}}{A list of length \code{nrow(dat.cluster)}. Each element is
#'   a \eqn{2 \times (N_i + 1)} matrix of binary potential outcomes for one
#'   individual. Rows correspond to \code{A_i = 0} and \code{A_i = 1}, and
#'   columns correspond to the number of treated neighbors.}
#'   \item{\code{dat.cluster}}{A data frame containing the clustered data
#'   together with:
#'   \describe{
#'     \item{\code{L_i}}{generated covariate, if needed}
#'     \item{\code{Ai}}{simulated individual treatment}
#'     \item{\code{treated_neighbors}}{number of treated neighbors in the cluster}
#'     \item{\code{Abar_Ni}}{proportion of treated neighbors}
#'     \item{\code{Y_obs}}{observed outcome under the realized treatment and neighbor exposure}
#'   }}
#' }
#'
#' @examples
#' dat.cluster <- build_cluster(K = 1000, lambda = 10)
#'
#' dat.cluster$X2 <- rbinom(nrow(dat.cluster), 1, 0.5)
#' dat.cluster$X3 <- runif(nrow(dat.cluster))
#' dat.cluster$X1 <- rbinom(nrow(dat.cluster), 1, 0.5)
#'
#' sim_out <- simulate_cluster_outcome(
#'   dat.cluster = dat.cluster,
#'   outcome_covars = c("L_i", "X1", "X2", "X3"),
#'   outcome_coefs = c(0.25, 0.6, -0.4, 0.2),
#'   treatment_covars = c("L_i", "X2"),
#'   treatment_coefs = c(1, 0.8)
#' )
#'
#' head(sim_out$PO, 10)
#' head(sim_out$dat.cluster)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate ungroup
#' @export

simulate_cluster_outcome <- function(dat.cluster,
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

  if ("L_i" %in% c(outcome_covars, treatment_covars) &&
      !("L_i" %in% names(dat.cluster))) {
    dat.cluster$L_i <- rbinom(nrow(dat.cluster), 1, p_L)
  }

  data <- dat.cluster
  n <- nrow(data)
  Ni <- data$interference_set_size

  if (length(outcome_covars) != length(outcome_coefs)) {
    stop("length(outcome_covars) must equal length(outcome_coefs)")
  }
  if (length(treatment_covars) != length(treatment_coefs)) {
    stop("length(treatment_covars) must equal length(treatment_coefs)")
  }

  if (!all(outcome_covars %in% names(data))) {
    stop("Some outcome_covars are not in dat.cluster")
  }
  if (!all(treatment_covars %in% names(data))) {
    stop("Some treatment_covars are not in dat.cluster")
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
  dat.cluster$Ai <- rbinom(nrow(dat.cluster), size = 1, prob = as.numeric(pA))

  dat.cluster <- dat.cluster %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::mutate(
      treated_neighbors = sum(Ai) - Ai,
      Abar_Ni = ifelse(interference_set_size == 0,
                       0,
                       treated_neighbors / interference_set_size)
    ) %>%
    dplyr::ungroup()

  dat.cluster$Y_obs <- vapply(
    seq_len(nrow(dat.cluster)),
    function(i) {
      PO_i <- PO[[i]]
      PO_i[
        dat.cluster$Ai[i] + 1,
        dat.cluster$treated_neighbors[i] + 1
      ]
    },
    numeric(1)
  )

  return(list(
    PO = PO,
    dat.cluster = dat.cluster
  ))
}

