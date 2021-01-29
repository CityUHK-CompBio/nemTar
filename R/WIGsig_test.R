#'  The statistical significance test based on the bootstrap of WIGs
#' @param sample_WIG0 the background distribution of the WIG from the sampling of WIG
#' @param WIG    the pathway-specific WIGcalculated using \code{compute_WIG}
#' @param sampling_times  a numeric number indicating the times of sampling
#' @return the statistical significance(BH-adjusted p-value)
#' @export
WIGsig_test<-function(sample_WIG0,WIG,sampling_times) {
  back_WIG<-sample_WIG0$back_WIGS
  TTT<-sample_WIG0$back_prob
  sig<-c()
  # EMT_WIGS<-test_WIGS$WIGS[c(2,6,7,8,9,11,12)]
  EMT_WIG<-WIG$WIG
  for (j in 1:length(EMT_WIG)){
    sig[j]<-length(which(back_WIG[,j]>=EMT_WIG[j]))/sampling_times
    names(sig)[j]<-names(EMT_WIG)[j]
  }
  p.adj<-stats::p.adjust(sig,method="BH")
  return(p.adj)
}
