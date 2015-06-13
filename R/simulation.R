#' Check matrix is symetric and positive semidefinate.
checkMatrix <- function(mat){
  if(any(eigen(mat)$values <= 0)) stop("Matrix is not positive semi-definate.")
}

#' Check ranges for eigen values are appropriate.
checkRange <- function(rng){
  if(any(rng <= 0)) stop("Range is not appropriate.")
}

#' Simulation of Multivariate Normal dataset
#' 
#' Creates a Multivariate Normal Dataset via random positive definate matrix generation.
#' 
#' @param num.vars The number of variables.
#' @param num.obs The number of observations.
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#' @return A data frame of multivariate normal observations.
#' @export
simMVNData <- function(num.vars, num.obs){ 
  cov.mat <- clusterGeneration::genPositiveDefMat(num.vars, covMethod=c("unifcorrmat"))$Sigma
  mvrnorm(num.obs, rep(0, num.vars), Sigma=cov.mat) %>%
    as.data.frame 
}

#' Random Graph Simulation
#' 
#' This is a wrapper for the bnlearn packages *random.graph* which gives a 
#' single graph with node names as numbers prefixed with "V", to match default
#' naming of variables in data.frames.
#' 
#' @param num.vars Number of variables
#' @param method A character string, possible values are \code{ordered}, 
#' \code{ic-dag}, and \code{melancon}. See \code{?random.graph} in \pkg{bnlearn}
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#' @return an object of class \code{bn}.
#' @export
randomNet <- function(num.vars, method){
  var.names <- paste("V", 1:num.vars, sep="")
  net <- num.vars %>%
{paste("V", 1:., sep="")} %>%
  random.graph(num = 1, method = method) 
}

#' Simulate a Sparce Precision Matrix
#' 
#' @param basis A binary matrix where 1's correspond to non-zero elements 
#' in the precision matrix. This is the output of the \code{getPrecisionBasis} function.
#' @param rng The desired range of the eigen values for the matrix.
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#' @return A positive definate matrix of the same dimension as the basis argument. 
#' @export
simSparsePrecision <- function(basis, rng){
  checkRange(rng)
  if(basis %>%
       as.numeric %>%
{!all(. %in% c(0, 1))}) stop("Basis elements must be 0 or 1.")
clusterGeneration::genPositiveDefMat(nrow(basis), 
                  covMethod="eigen",
                  lambdaLow = rng[1],
                  ratioLambda = rng[2] / rng[1]) %$% 
  Sigma %>%
  `*`(basis)
}

#' Get the precision matrix structure for a Gaussian joint distribution based on a Bayesian 
#' network structure
#'
#' Converts a Bayesian network to a binary matrix where 1's correspond to non-zero elements 
#' in the precision matrix. 
#' 
#' @param net A object of class \code{bn}.
#' 
#' @return A symetric matrix.
#' @export
getPrecisionBasis <- function(net){
  net %>% 
    moral %>%
    amat %>%
    `diag<-`(1)
}

#' Simulate a covariance matrix based on a Bayesian network structure
#' 
#' Given a Bayesian network structure, simulates a covariance matrix consistent 
#' with the network structure, based on Gaussian assumptions.
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#' @param net A Bayesian network structure, an object of class bn.
#' @param rng The desired range of the eigen values in the covariance matrix.
#' @return A simulated the covariance matrix consistent with the net argument. Simulated 
#' under the Gaussian assumption.
#' @export
simNetCovariance <- function(net, rng){
  checkRange(rng)
  getPrecisionBasis(net) %>% 
    simSparsePrecision(rng[c(2, 1)]^-1) %>%
    solve
}
#' Convert to Binary Matrix
#' 
#' Converts a matrix to a binary matrix, where 0 values remain 0 values, and non-zero values 
#' become 1.
convert2BinaryMat <- function(mat){
  if(length(dimnames(mat)) == 0){
    var.names <- paste("V", 1:nrow(mat), sep = "")
    dimnames(mat) <- list(var.names, var.names)
  }
  (mat != 0) * 1
}

#' Convert Adjacency matrix to Moral (Undirected Conditional Dependence) Network
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#' @param adj.mat An adjacency matrix
#' @return An object of class bn where all edges are undirected
adjMat2Moral <- function(adj.mat){
  base.net <- adj.mat %>% colnames %>% empty.graph 
  amat(base.net, ignore.cycles=T) <- adj.mat
  base.net
}

#' Simulate a Gaussian Bayesian network structure based on a covariance matrix.
#' 
#' This function attempts to find a directed acyclic structure consistent with the
#' conditional independencies given by the covariance matrix under Gaussian assumptions.
#' If such a consistent extention cannot be found, then NULL is returned.
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'   \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#' @param cov.mat A covariance matrix (positive definate matrix).
#' @return An objective of class bn.  If no consistent extention is found, NULL is returned.
#' @export
simDAGFromCovMat <- function(cov.mat){
  checkMatrix(cov.mat)
  tryCatch(cov.mat %>%
             solve %>%
             round(4) %>% 
             convert2BinaryMat %>% 
             `diag<-`(0) %>%
             adjMat2Moral %>%
             cextend ,
           error = function(e) NULL)
}

#' Simulate a Bayesian network from a covariance matrix
#' 
#' Given a covariance matrix, simulates a fully parameterized Guassian Bayesian network 
#' consistent with the covariance matrix.
#' 
#'  @param cov.mat A positive definate matrix. 
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#'  @return An object of class \code{bn.fit}. If no consistend extention is found, 
#'  NULL is returned.
#'  @export
simBNFromCovMat <- function(cov.mat){   
  structure <- simDAGFromCovMat(cov.mat) #Get a DAG from the covariance matrix
  if(is.null(structure)) return(NULL)
  sim.data <- mvrnorm(1000, rep(0, ncol(cov.mat)), Sigma=cov.mat) %>%
    as.data.frame
  bn.fit(structure, sim.data, method = "mle")  
}

#' Simulate a fully specified Gaussian Bayesian network
#' 
#' Simulates a Gaussian Bayesian network by simulating a network structure and 
#' a random covariance matrix consistent with the I-map of the network structure, 
#' and then calculating network parameters based on the covariance matrix.
#' 
#' @param num.vars Number of variables/nodes
#' @param rng The desired range of the eigen values in the covariance matrix.
#' @param method A character string, possible values are \code{ordered}, 
#' \code{ic-dag}, and \code{melancon}. See \code{?random.graph} in \pkg{bnlearn}
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#' @return An object of class \code{bn.fit} 
#' @export
simGaussianNet <- function(num.vars, rng = c(.5, 1.5), method = "melancon"){
  checkRange(rng)
  random.net <- randomNet(num.vars, method) 
  cov.mat <- random.net %>%  simNetCovariance(rng) 
  # simBNFromCovMat relies on cextend, which may not work for dense graphs.
  # In this case, simBNFromCovAndNet is used.
  sim.net <- simBNFromCovMat(cov.mat)
  if(is.null(sim.net)) {
    sim.net <- simBNFromCovAndDAG(random.net, cov.mat)
  }
  sim.net
}
#' Simulate a Bayesian network from a DAG and covariance matrix
#' 
#' Given a covariance matrix and a directed acyclic graph structure, 
#' simulates a fully parameterized Guassian Bayesian network consistent 
#' with the covariance matrix.
#' 
#'  @seealso \code{\link{simGaussianNet}}, \code{\link{simBNFromCovMat}}, 
#'  \code{\link{simDAGFromCovMat}}, \code{\link{simNetCovariance}}, \code{\link{simMVNData}},
#'  \code{\link{randomNet}}, \code{\link{simSparsePrecision}}, \code{\link{simDAGFromCovMat}}
#'  
#'  @param dag A object of class bn, with no undirected edges.
#'  @param cov.mat A positive definate matrix. 
#'  @return An object of class \code{bn.fit}
#'  @export
simBNFromCovAndDAG <- function(dag, cov.mat){
  if(nrow(undirected.arcs(dag)) > 0) stop("Input DAG cannot have undirected arcs.")
  sim.data <- mvrnorm(1000, rep(0, ncol(cov.mat)), Sigma=cov.mat) %>%
    as.data.frame
  bn.fit(dag, sim.data, method="mle")
}
