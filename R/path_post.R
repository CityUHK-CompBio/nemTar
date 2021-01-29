#' To dissect the influenced E-genes and the probalility of influence under the assumption
#' of nested effects for a specific pathway
#'
#' @param nemTar_rslt the output of the signaling network inferred by nemTar
#'
#' @param path_list  a character of the the signature genes of an pathway

#' @param Sgenes  a character of the S-genes within the signaling network
#'
#' @import dplyr
#' @export
#############################################################################
path_post<-function(nemTar_rslt,path_list,Sgenes) {
  Sig_stat<-methods::as(nem::transitive.closure(nemTar_rslt$graph),"matrix")
  sig_list<-path_list
  # pathway_related affected genes extraction
  path_related<-list()
  Sig_affected<-list()
  path_related<-lapply(nemTar_rslt$mappos,function(s) intersect(s,sig_list))
  Sig_affected<-apply(Sig_stat,1,function(s) which(s==1))
  path_affected<-list()
  # for (i in 1:length(Sig_affected)) {
  #   path_affected[[i]]<-sapply(path_related[unlist(Sig_affected[[i]])],function(s) unique(unlist(s)))
  # }
  # path_affected <- lapply(1:length(Sig_affected), function(x) {
  #   unique(unlist(path_related[unlist(Sig_affected[[x]])]))
  #
  # })
  path_affected <- lapply(1:length(Sig_affected), function(x) {
    unlist(path_related[unlist(Sig_affected[[x]])])%>%unique
  })
  names(path_affected)<-rownames(Sig_stat)

  posprob_affected<-list()
  length(posprob_affected)<-length(Sgenes)
  names(posprob_affected)<-names(path_affected)
  leng<-c()
  for (m in 1:length(path_affected)){
    leng[m]<-length(path_affected[[m]])
  }
  num<-which(leng!=0)
  for (i in num){
    for (j in 1:length(path_affected[[i]])){
      if (i %in% grep(path_affected[[i]][j],path_related)){
        posprob_affected[[i]][j]<-nemTar_rslt$pos[match(path_affected[[i]][j],rownames(nemTar_rslt$pos)),i]
      }
      else{
        posprob_affected[[i]][j]<-max(nemTar_rslt$pos[match(path_affected[[i]][j],rownames(nemTar_rslt$pos))
                                           ,grep(path_affected[[i]][j],path_related)])
      }
    }
  }
  return(list(Sig_affected=Sig_affected,path_affected=path_affected,post_affected=posprob_affected))
}



