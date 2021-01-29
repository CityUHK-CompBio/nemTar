#' Calculating the weigheted information gain(WIG) of combinational perturbation
#' of 2 S-genes
#'
#' @param path_post the influenced E-genes and the probalility of influence
#'   under the assumption of nested effects
#' @param num_sgene      a numeric of the number of S-genes
#'
#' @import dplyr
#' @export
WIG_double<-function(path_post,num_sgene){
  ######## 2 S-genes pertubation
  multi_s0<-path_post$path_affected
  multi_posprob0<-path_post$post_affected
  for (i in 1:length(multi_s0)){
    if (length(multi_s0[[i]])!=0)
      names(multi_posprob0[[i]])<-multi_s0[[i]]
  }
  # names(multi_posprob0[[1]])
  overlap2_Sgene<-apply(t(combn(c(1:length(multi_s0)),2)),1,function(x)
    union(multi_s0[[x[1]]],multi_s0[[x[2]]]))
  leng_overlap2_Sgene<-dim(t(combn(c(1:length(multi_s0)),2)))[1]
  for (k in 1:leng_overlap2_Sgene){
    names(overlap2_Sgene)[k]<-paste0(names(multi_s0)[t(combn(c(1:length(multi_s0)),2))[k,1]],"/",
                                     names(multi_s0)[t(combn(c(1:length(multi_s0)),2))[k,2]])
  }
  # overlap2_egene_posprob<-apply(t(combn(c(1:length(multi_s0)),2)),1,function(x)
  #   union(multi_posprob0[[x[1]]],multi_posprob0[[x[2]]]))
  # names(overlap2_egene_posprob)<-names(overlap2_sgene)
  overlap2_Egene_posprob<-apply(t(combn(c(1:length(multi_s0)),2)),1,function(x)
    ## compare whether the intersection exist
    if ((intersect(multi_s0[[x[1]]],multi_s0[[x[2]]]) %>% length)==0){
      c(multi_posprob0[[x[1]]],multi_posprob0[[x[2]]])
    }
    else {
      overlapped<-intersect(multi_s0[[x[1]]],multi_s0[[x[2]]])
      ##find the overlapped E-genes in x[[1]] and x[[2]]
      num_set1<-match(overlapped,multi_s0[[x[1]]])
      num_set2<-match(overlapped,multi_s0[[x[2]]])
      #3 only 1 and 2
      only1<-setdiff(multi_s0[[x[1]]],multi_s0[[x[2]]])
      only2<-setdiff(multi_s0[[x[2]]],multi_s0[[x[1]]])
      com_diff<-multi_posprob0[[x[1]]][num_set1]-multi_posprob0[[x[2]]][num_set2]
      c(multi_posprob0[[x[1]]][match(only1,names(multi_posprob0[[x[1]]]))],multi_posprob0[[x[2]]][match(only2,names(multi_posprob0[[x[2]]]))],
        multi_posprob0[[x[1]]][num_set1][which(com_diff>=0)],multi_posprob0[[x[2]]][num_set2][which(com_diff<0)])
    }
  )
  names(overlap2_Egene_posprob)<-names(overlap2_Sgene)
  WIG_2Sgene<-compute_WIG(overlap2_Egene_posprob,overlap2_Sgene,num_sgene)
  #sort(WIGS_2Sgene$WIG,decreasing = T)
  return(WIG_2Sgene=WIG_2Sgene)
}
