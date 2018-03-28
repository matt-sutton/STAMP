#' Simulated Dataset SC11
#'
#' Example dataset from the simulation study in the paper.
#'
#' @docType data
#'
#' @usage data(Example_SC11)
#'
#' @keywords datasets
#'
#' @references Andriy Derkach, Ruth M. Pfeiffer (2018)
#' (\href{https://arxiv.org/abs/1709.06702}{Arxiv})
#'
#' @examples
#' data(Example_SC11)
#' # run mixure model with quadratic statistic
#' \donttest{mixtureQuad(1000,Example_SC11$ListBeta,Example_SC11$ListSE,Example_SC11$ListCor,Example_SC11$MafList) }
"Example_SC11"
