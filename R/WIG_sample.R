#' Sampling of the WIG to measure the statistical significance
#'
#' @param path_post the influenced E-genes and the probalility of influence
#'   under the assumption of nested effects
#' @param nem_rslt    the output of the signaling network inferred by nemTar
#' @param num_sgene      a numeric of the number of S-genes
#' @param sampling_times a numeric number indicating the times of sampling
#' @param path_list a character of the the signature genes of an pathway
#'
#' @import dplyr
#' @export

WIG_sample<-function(path_post,nem_rslt,path_list,num_sgene,sampling_times){
  mappos_E <- sapply(1:length(path_post$Sig_affected), function(x) {
    unlist(nem_rslt$mappos[path_post$Sig_affected[[x]]])%>%unique
  })
  names(mappos_E)<-names(path_post$Sig_affected)

  pos_E<-sapply(1:length(mappos_E),function(i) nem_rslt$pos[match(mappos_E[[i]],rownames(nem_rslt$pos)),])
  names(pos_E)<-names(path_post$Sig_affected)
  prob_E0<-list()
  for (k in 1:length(pos_E)){
    prob_E0[[k]]<-apply(pos_E[[k]],1,function(s) s[k])
    path_id<-match(intersect(rownames(pos_E[[k]]),path_list),rownames(pos_E[[k]]))
    if (length(path_id)>1)
      prob_E0[[k]][path_id]<-apply(pos_E[[k]][path_id,],1,max)
    else if(length(path_id)==1)
      prob_E0[[k]][path_id]<-max(pos_E[[k]][path_id,])
  }
  names(prob_E0)<-names(mappos_E)
  ##########################################################################################################
  ##########################################################################################################
  num_choose<-c()
  for (k in 1:length(path_post$path_affected)){
    if (lengths(path_post$path_affected)[k]!=0){
      num_choose[k]<-choose(lengths(mappos_E)[k],lengths(path_post$path_affected)[k])
    }
    else {
      num_choose[k]<-0
    }
  }
  # remove the S-genes that have the # of maximum choose of combination of Egenes less than given sample times and the "null"S-gene
  sample_mappos<-mappos_E[which(lengths(path_post$path_affected)!=0)]
  sample_posprob<-prob_E0[which(lengths(path_post$path_affected)!=0)]
  sample_Egene_num<-lengths(path_post$path_affected)[which(lengths(path_post$path_affected)!=0)]
  back_prob<-vector("list",sampling_times)
  back_prob<-lapply(1:sampling_times, function(i){
    # res1<-lapply(names(sample_mappos), function(j){
    #   sample_posprob[[j]][sample(sample_mappos[[j]], sample_Egene_num[j],replace=T)]
    # })
    res1<-list()
    for (j in 1:length(sample_mappos)){
      res1[[j]]<-sample_posprob[[j]][match(sample(sample_mappos[[j]], sample_Egene_num[j],replace=T),
                                           sample_mappos[[j]])]
    }
    return(res1)
  })
  ##### Compute the background distribution of WIGS
  p0<-1/(num_sgene+1) #(+1 for a "null" S-gene)
  # back_WIGS0<-Map(function(x,y) back_prob[[x]][[y]]*log(back_prob[[x]][[y]]/p0),
  #                    x=1:sampling_times,y=1:length(sample_mappos))
  # back_WIGS<-mapply(function(x,y) sum(back_WIGS0[[x]][[y]]),
  #                   x=1:sampling_times,y=1:length(sample_mappos)) # some strtegies need trial
  back_WIGS0<-vector("list",sampling_times)
  back_WIGS<-matrix(NA,nrow=sampling_times,ncol=length(sample_mappos))
  colnames(back_WIGS)<-names(sample_mappos)
  # for (i in 1:sampling_times){
  #   for (j in 1:length(sample_mappos)){
  #     back_WIGS0[[i]][[j]]<-sapply(back_prob[[i]][[j]],function(s) s*log(s/p0))
  #   }
  # }
  back_WIGS0<-lapply(1:sampling_times, function(i){
    res2<-lapply(back_prob[[i]], function(x){
      x*log(x/p0)
    })
    return(res2)
  })
  for (i in 1:sampling_times){
    back_WIGS[i,]<-sapply(back_WIGS0[[i]], function(t) sum(t)) # some strtegies need trial
  }
  return(list(num_choose=num_choose,back_prob=back_prob,back_WIGS0=back_WIGS0,back_WIGS=back_WIGS))
}
