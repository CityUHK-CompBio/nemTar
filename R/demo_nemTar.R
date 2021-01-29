


# demo_nemTar<-function(method,Sgenes,S_mat,E_mat,para,sig_genes,sample_times){
#   if (method=="gre"){
#     control<-nem::set.default.parameters(Sgenes,type="mLL",para=para)
#     nemTar_rslt<-nemTar_greedy(E_mat,Sgenes,S_mat,control=control)
#   }
#   else if (method == "tri"){
#   control<-nem::set.default.parameters(Sgenes,type="mLL",para=para)
#   nemTar_rslt<-nemTar_triples(E_mat,Sgenes,S_mat,control=control)
#   }
#   else
#  {cat("Sorry,nemTar cannot support other unspecified inference methods!")}
#   if (method=="gre"|method == "tri"){
#    eff_path<-path_post(nemTar_rslt,sig_genes)
#    WIG_rslt<-compute_WIG(eff_path$posprob_affected,eff_path$path_affected,length(Sgenes))
#    sample_WIG0<-sample_WIG(nemTar_rslt,eff_path$path_affected,length(Sgenes),sample_times)
#    back_WIGS<-sample_WIG0$back_WIGS
#    sig<-c()
#    path_WIG<-WIG_rslt$WIG
#    for (j in 1:length(path_WIG)){
#      sig[j]<-length(which(back_WIGS[,j]>=path_WIG[j]))/sample_times
#    }
#    p_adj<-p.adjust(sig,method="BH")
#    names(p_adj)<-names(path_WIG)
#    WIG_output<-rbind(path_WIG, WIG_output)
#   }
#   return(list(nemTar_rslt=nemTar_rslt,path_WIG=path_WIG,WIG_output=WIG_output))
# }
