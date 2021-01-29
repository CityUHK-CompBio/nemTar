#' Employing the permuation test
#' @param path_post the influenced E-genes and the probalility of influence
#'   under the assumption of nested effects
#' @param nem_rslt   the output of the signaling network inferred by nemTar
#' @param path_list  a character of the the signature genes of an pathway
#' @param sampling_times  a numeric number indicating the times of sampling
#' @import dplyr
#'
permu_test<-function(path_post,nem_rslt,path_list,sampling_times) {
  mappos_E <- sapply(1:length(path_post$Sig_affected), function(x) {
    unlist(nem_rslt$mappos[path_post$Sig_affected[[x]]])%>%unique
  })
  names(mappos_E)<-names(path_post$Sig_affected)

  EMT_sig<-unique(unlist(path_post$path_affected))
  draw_E<-mappos_E[which(lengths(path_post$path_affected)!=0)]
  E_all<-unique(unlist(nem_rslt$mappos))

  draw_EMT<-path_post$path_affected[which(lengths(path_post$path_affected)!=0)]
  #sampling_times<-10000
  back_sam<-vector("list",sampling_times)
  # back_sam<-lapply(1:sampling_times, function(i){
  #   res1<-lapply(draw_EMT, function(j){
  #     EMT_sig[sample(c(1:length(EMT_sig)),lengths(draw_EMT)[j],replace=T)]
  #   })
  #   return(res1)
  # })
  for (i in 1:sampling_times){
    for (j in 1:length(draw_EMT)){
      back_sam[[i]][[j]]<-E_all[sample(c(1:length(E_all)),lengths(draw_E)[j],replace=T)]
    }
  }
  retri_mat<-matrix(NA,nrow=sampling_times,ncol=length(draw_EMT),dimnames = list(c(1:sampling_times),names(draw_EMT)))
  for (i in 1:sampling_times){
    for (j in 1:length(draw_EMT)){
      retri_mat[i,j]<-length(intersect(back_sam[[i]][[j]],path_list))
    }
  }
  sig<-c()
  for (j in 1:length(draw_EMT)){
    sig[j]<-length(which(retri_mat[,j]>=lengths(draw_EMT)[j]))/sampling_times
    names(sig)[j]<-names(draw_EMT)[j]
  }
  p.adj<-stats::p.adjust(sig,method="BH")
  #rslt<-data.frame(S_genes=names(p.adj),adjusted_P=p.adj)
  return(p.adj)
}
