#' Calculating the weigheted information gain(WIG) of each S-genes
#'
#' @param post_prob a list of the probablitiy that an E-gene is affected by an
#'   S-gene with respect to a specific pathway
#' @param path_affected  a list of the E-genes that are affected by the S-genes
#'   and are realted to a specific pathway
#' @param num_sgene  a numeric of the number of S-genes
#'
#'
#' @return a list contains the probability of the uniform distribution before
#'   the nertwork inference, and the WIG profiles of every S-gene,as well as the
#'   WIGs of S-genes
#' @export

compute_WIG<-function(post_prob,path_affected,num_sgene){
  p0<-1/(num_sgene+1) #(+1 for a "null" S-gene)
  len<-c()
  for (m in 1:length(path_affected)){
    len[m]<-length(path_affected[[m]])
  }
  WIG0<-sapply(post_prob[which(len!=0)],function(s) s*log(s/p0))
  #WIG<-sapply(WIG0, function(t) 0.5*max(t)+0.5*median(t))
  WIG<-sapply(WIG0, function(t) sum(t)) # other strtegies could also be tested
  return(list(p0=p0,WIG0=WIG0,WIG=WIG))
}
