#' Compute the concetration matrix via IPS.
#'
#' @param p A correlation between 0 and 1.
#' @param A An adjacency matrix.
#' @param eps (optional) Convergence tolerance.
#' @return The concentration matrix \code{Q}.
#' @examples
#' source("http://www.biostat.umn.edu/~hodges/RPLMBook/Datasets/10_Slovenian_stomach_cancer/Slovenia_stomach_cancer_data.txt")
#' A       <- -Q
#' diag(A) <- 0
#' Q       <- IPS(0.5, A)
#'
#' @seealso \code{\link{clique_part}}
#' @export
IPS <- function(p, A, eps = 1e-8){
  n <- nrow(A)
  diag(A) <- 0
  T <- diag(n) + p*A

  W <- clique_part(A)
  Q <- IPS2(rho=p, n=n, T=T, W=W, maxit=100, eps=eps)$Q

  return(Q)
}



#' Generate a clique partition of the graph G defined by the adjacency matrix A.
#'
#' @param A An adjacency matrix.
#' @return A clique partition in a list.
#' @examples
#' clique_part(A)
clique_part <- function(A){
  G      <- igraph::graph_from_adjacency_matrix(A)
  Clist  <- suppressWarnings(igraph::max_cliques(G))
  Clist  <- lapply(Clist, function(x) sort(as.numeric(x)))

  V      <- sum(lower.tri(A)*A)
  NN     <- length(Clist)
  Nc     <- NN/(1:NN)
  NcTemp <- Nc[which(Nc %% 1 == 0)]
  M      <- NcTemp[-c(1, length(NcTemp))]

  Mcomp  <- sapply(M, clique_select,Clist=Clist, V=V, NN=NN)

  nc     <- M[which(Mcomp==min(Mcomp))]
  nC     <- length(Clist)/nc
  W      <- lapply(1:nc, function(i) Clist[(i*nC-nC+1):(i*nC)])

  #Need to adjust index since C++ index begins with zero
  W1  <- lapply(1:nc, function(i) lapply(W[[i]], function(x) x-1))
  W1
}


#' Compute the computational complexity associated with a clique partition.
#'
#' @param m Size of the parition.
#' @param Clist A list of cliques.
#' @param V Number of vertices in the graph.
#' @param NN Number of cliques in the graph.
#' @return The computational complexity associated with a clique partition.
#' @examples
#' #Internal function only
clique_select <- function(m, Clist, V, NN){
  nC  <- length(Clist)/m
  Wt  <- lapply(1:m, function(i) Clist[(i*nC-nC+1):(i*nC)])
  Umstar <- max(sapply(1:m, function(j) sum(sapply(Wt[[j]], length))))
  log( m*V^3 + NN * Umstar^3 )
}
