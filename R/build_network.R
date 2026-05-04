#' Build a Network from Ego-Alter Edge Lists
#'
#' This function constructs an undirected network using ego and alter IDs.
#' It removes isolated vertices, simplifies the graph by removing duplicate
#' edges and self-loops, and returns the resulting graph along with node
#' degree information.
#'
#' @param x1 A vector containing ego node IDs.
#' @param x2 A vector containing alter node IDs.
#'
#' @return A list containing:
#' \describe{
#'   \item{net0}{An igraph network object.}
#'   \item{n}{Number of nodes in the network.}
#'   \item{data}{A data frame with node IDs and their degrees.}
#' }
#'
#' @examples
#' ego <- c(1,1,2,3,4,5)
#' alter <- c(2,3,3,4,5,1)
#' res <- build_network(ego, alter)
#' res
#'
#' @export




#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'



build_network <- function(x1, x2) {

  edges <- data.frame(
    EGO_ID   = as.character(x1),
    ALTER_ID = as.character(x2)
  )

  nodes <- data.frame(
    name = unique(c(edges$EGO_ID, edges$ALTER_ID))
  )

  net0 <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )

  net0 <- igraph::delete_vertices(net0, which(igraph::degree(net0) == 0))
  net0 <- igraph::simplify(net0, remove.multiple = TRUE, remove.loops = TRUE)

  n <- igraph::vcount(net0)

  data <- data.frame(
    subject_id = as.integer(igraph::V(net0)$name),
    id = seq_len(n)
  )

  data$na <- igraph::degree(net0)

  return(list(
    net0 = net0,
    n = n,
    data = data
  ))
}
