#' Bootstrap Learning Performed on a Gaussian Bayesian Network.
#'
#'  The 15 node network structure was simulated with Melancon's and Philippe's
#'  Uniform Random Acyclic Digraphs algorithm.  This produces a highly connected DAG.
#'  See \code{?random.graph} for more information.
#'
#' \itemize{
#'   \item truth. An object of \code{bn.fit}, a simulated structure and parameter set.
#'   \item boot.  An object of \code{bn.strength} and \code{data.frame}, results of
#'   model averaging applied to data simulated from truth.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name melancon_boot
#' @usage data(melancon_boot)
#' @format A list with 2 items
NULL

#' Bootstrap Learning Performed on a Gaussian Bayesian Network.
#'
#'  The 40 node network structure was simulated with full ordering based generation.
#'  This produces a sparsely connected DAG.  See \code{?random.graph} for more information.
#'
#' \itemize{
#'   \item truth. An object of \code{bn.fit}, a simulated structure and parameter set.
#'   \item boot.  An object of \code{bn.strength} and \code{data.frame}, results of
#'   model averaging applied to data simulated from truth.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ordered_boot
#' @usage data(ordered_boot)
#' @format A list with 2 items
NULL
