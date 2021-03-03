#' Generation of the data for sampling in the simulation of the cancer multi-omics samples
#'
#' @param Phi adjacency matrix
#' @param S_obs  binary matrix inddicating the observation of S-genes' states
#' @param p  number of samples(patients) to sample
#' @param m  number of E-genes to sample
#' @param type   an integer,represents the number of row within the raw image
#' @param typeI.err   simulated type I error for binary data
#' @param typeII.err  simulated type II error for binary data

#' @seealso \code{\link[nem]{sampleData}}
sampleData2 = function(Phi,S_obs,p,m,type="binary", typeI.err=typeI.err, typeII.err=typeII.err){
  Sgenes = colnames(Phi)
  n = length(Sgenes)
  set.seed(1234)
  epos = sample(1:n,m,replace=TRUE,prob=prob)
  Theta = matrix(0, nrow=length(Sgenes), ncol=m)
  for(i in 1:ncol(Theta))
    Theta[epos[i],i] = 1
  Phi = transitive.closure(Phi, mat=TRUE, loops=TRUE)
  M = t((S_obs %*% Phi%*%Theta > 0)*1)
  if(type == "binary"){
    D = matrix(0, ncol=p, nrow=m)
    k2 = 1
    for(i in 1:p){
      D[M[,i] == 1, k2] =  matrix(sample(c(0,1),sum(M[,i] == 1),replace=TRUE,prob=c(typeII.err,1-typeII.err)),ncol=1)# effected genes => aus H1 ziehen
      D[M[,i] == 0, k2] =  matrix(sample(c(0,1),sum(M[,i] == 0),replace=TRUE,prob=c(1-typeI.err,typeI.err)), ncol=1) # ... and not effected ones => aus H0 ziehen
      k2 = k2 + 1
    }
    # if(uninformative > 0)
    #   D = rbind(D, matrix(sample(c(0,1),n*uninformative,replace=TRUE,prob=c(1-typeI.err,typeI.err)),
    #                       nrow=uninformative, ncol=n)) # ... and not effected ones => aus H0 ziehen)
  }
  # else if(type %in% c("density")){
  #   lambda = cbind(1-rowSums(lambda), lambda)
  #   palt = sapply(1:n, function(i) bum.ralt(m, c(alpha[i], beta[i]), lambda[i,]))
  #   p0 = matrix(runif((m+uninformative)*n), ncol=n)
  #   P = M*palt + (1-M)*p0[1:m,]
  #   if(uninformative > 0)
  #     P = rbind(P, p0[(m+1):nrow(p0),])
  #   D = sapply(1:n, function(i) bum.dalt(P[,i], c(alpha[i], beta[i]), lambda[i,]))
  #   D = log(D)
  # }
  # else if(type %in% c("lodds")){
  #   palt = sapply(1:n, function(i) rnorm(m, mean=meansH1[i], sd=sdsH1[i]))
  #   p0 = sapply(1:n, function(i) rnorm(m+uninformative, mean=meansH0[i], sd=sdsH0[i]))
  #   D = M*palt + (1-M)*p0[1:m,]
  #   if(uninformative > 0)
  #     D = rbind(D, p0[(m+1):nrow(p0),])
  # }
  else
    stop(paste("unknown type", type, "\n"))
  if(type == "binary")
    colnames(D) = rownames(S_obs)
  else
    colnames(D) = Sgenes
  list(D=D, epos=Sgenes[epos])
}
