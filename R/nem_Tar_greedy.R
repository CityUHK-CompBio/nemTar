#' Inferring the signaling network from (epi)genetic aberrations and effect
#' reporter gene profiles based on greedy hill climbing methods
#'
#' @param D an E-gene observation matrix with the dimension of p(number of
#'   samples) × n(number of effect reporter genes)
#' @param Sgenes  a character vector of S-genes
#' @param S_pattern  a binary matrix showing the state of whether signaling
#'   components(S-genes)are perturbed, with the dimension of m × m
#' @param initial the initial state of the S-gene graph
#' @param control a list contains necessary parameters for implementing the
#'   methods,the details could be referred to \code{\link[nem]{nem}}
#' @param verbose a logic to whether present the execution process
#' @return a list contains the network inference results, please also refer to
#'   \code{\link[nem]{nem}}
#' @import nem
#' @export
#' @seealso \code{\link{nem_Tar_triples}} \code{\link[nem]{nem}}

nem_Tar_greedy <- function(D,Sgenes,S_pattern,initial=NULL,control, verbose=TRUE){
  # if(is(D, "list"))
  #   Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D[[1]]))]),"time")
  # else
  #   Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
  n <- length(Sgenes)
  cat("Greedy hillclimber for",n,"S-genes (lambda =", control$lambda,")...\n\n")
  if(is.null(initial))
    Phi <- matrix(0,nrow=n,ncol=n)
  else
    Phi = initial
  diag(Phi) <- 1
  dimnames(Phi) <- list(Sgenes,Sgenes)
  sco0 <-search2(D,Sgenes,S_pattern,control,models=list(Phi))$mLL
  finished <- FALSE
  while(!finished){
    models <- list()
    # propose new edges
    models = get.insertions(Phi, control$trans.close)
    if(control$type %in% c("CONTmLLMAP", "CONTmLLRatio") & !control$trans.close){ # these graphs are NOT transitively closed necessarily
      models = c(models, get.deletions(Phi), get.reversions(Phi))
    }
    models <- unique(models)
    if(verbose)
      cat(length(models), " local models to test ...\n")
    if(length(models) > 0){
      sconew <-
        search2(D,Sgenes,S_pattern,control,models)
      if(max(sconew$mLL) > sco0){
        if(verbose)
          cat("--> Edge added, removed or reversed\n")
        sco0 <- max(sconew$mLL)
        Phi <- methods::as(sconew$graph,"matrix")
      }
      else # otherwise no improving edge could be inserted
        finished <- TRUE
    }else
      finished <- TRUE
  }
  if(control$backward.elimination & !control$trans.close){
    if(verbose)
      cat("Backward elimination step:\n\n")
    finished <- FALSE
    while(!finished){
      # delete edges
      idx = which(Phi - diag(n) == 1)
      if(length(idx) > 0){
        models <- list()
        for(i in 1:length(idx)){ # test all possible deletions
          Phinew = Phi
          Phinew[idx[i]] = 0
          models[[i]] <- Phinew
        }
        models <- unique(models)
        sconew <- nem(D, models=models, inference="search", control=control, verbose=verbose)
        if(max(sconew$mLL) > sco0){
          if(verbose)
            cat("--> Edge deleted\n")
          sco0 <- max(sconew$mLL)
          Phi <- methods::as(sconew$graph,"matrix")
        }
        else # otherwise no improving edge could be deleted
          finished <- TRUE
      }else
        finished <- TRUE
    }
  }
  ep <- search2(D,Sgenes,S_pattern,control,models=list(Phi))
  res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control,
              selected=ep$selected, LLperGene=ep$LLperGene[[1]], para=ep$para[[1]])  # output: data likelihood under given model!
  class(res) <- "nem.greedy"
  if(verbose)
    cat("log-likelihood of model = ",res$mLL,"\n")
  return(res)
}

## Performing the insertion of edges in greedy hill climbing methods
get.insertions = function(Phi, trans.close=TRUE){
  idx = which(Phi == 0)
  models = list()
  if(length(idx) > 0){
    for(i in 1:length(idx)){ # test all possible new edges
      Phinew = Phi
      Phinew[idx[i]] = 1
      if(trans.close)
        Phinew = transitive.closure(Phinew, mat=TRUE,loops=TRUE)
      models[[i]] <- Phinew
    }
  }
  models
}

## Performing the deletion of edges in greedy hill climbing methods
get.deletions = function(Phi){
  Phi = Phi - diag(ncol(Phi))
  idx = which(Phi == 1)
  models = list()
  if(length(idx) > 0){
    for(i in 1:length(idx)){ # test all possible edge deletions
      Phinew = Phi
      Phinew[idx[i]] = 0
      diag(Phinew) = 1
      models[[i]] <- Phinew
    }
  }
  models
}

## Performing the reversion of edge direction in greedy hill climbing methods
get.reversions = function(Phi){
  idx = which(Phi + t(Phi) == 1, arr.ind=TRUE)
  models = list()
  if(NROW(idx) > 0){
    for(i in 1:NROW(idx)){ # test all possible edge reversions
      Phinew = Phi
      Phinew[idx[i,1],idx[i,2]] = 0
      Phinew[idx[i,2],idx[i,1]] = 1
      diag(Phinew) = 1
      models[[i]] <- Phinew
    }
  }
  models
}
