## To search for the best model based on marginal log likelihood
search2<-function(D,Sgenes,S_pattern,control,models){
  D<-D
  Sgene_obs<-S_pattern
  sam<-rownames(Sgene_obs)
  # colnames(D)<-sam
  D1 = sapply(sam, function(s) rowSums(D[,colnames(D) == s,drop=FALSE]))
  D0 = sapply(sam, function(s) sum(colnames(D) == s)) - D1

  S_state<-list()
  S_state<-lapply(models, function(s) sign(Sgene_obs %*% s))

  results <- lapply(S_state,mLL2,D1,D0,control,verbose=TRUE)
  # summerize the inference results
  s = sapply(results, function(r) r$mLL)
  ep   <- lapply(results, function(r) r$pos)
  map = lapply(results, function(r) r$mappos)
  LLperGene = lapply(results, function(r) r$LLperGene)
  para = lapply(results, function(r) r$para)
  selected = map[[which.max(s)]]
  selected = unique(unlist(selected[Sgenes]))
  if(length(s) > 1){
    mLL.sorted = sort(s, decreasing=TRUE)
    cat("((Marginal) posterior likelihood difference of best vs. second best model for ", Sgenes, ":", mLL.sorted[1] - mLL.sorted[2],")\n")
  }
  # winning model
  winner <- models[[which.max(s)]]
  diag(winner) <- 0
  gR <- methods::new("graphAM",adjMat=winner[Sgenes,Sgenes],edgemode="directed")
  gR <- methods::as(gR,"graphNEL")
  res <- list(graph=gR, mLL=s, pos=ep, mappos=map, control=control, selected=selected, LLperGene=LLperGene, para=para)
  class(res) <- "score"
  return(res)
}
########################################################################################################
