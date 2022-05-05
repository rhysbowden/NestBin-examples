library("NestBin")
#' Simulate data from a trial with a stepped wedge design 
#'
#' Function to sample binary values from a longitudinal cluster randomized 
#'   design with given prevalences and a block-exchangeable correlation 
#'   structure.
#'
#' @param design_matrix is 0-1 matrix with one row per sequence and one column per time 
#'   period. A value of 1 indicates an additive treatment effect of treatment_effect.
#' @param time_effect_matrix is the same dimensions as design_matrix, but it represents 
#'   the prevalence in each sequence-period cell, excluding the treatment effect.
#' @param icc is the within-period (cluster-period cell level) correlation.
#' @param cac is the cluster auto-correlation
#' @param n is number of observations per cluster-period
#' @param C is number of clusters per sequence
#' @return a data frame with the following columns:
#'   outcome: the 0-1 outcomes
#'   C_labels,T_labels,CT_labels: integers representing the cluster, period or 
#'       cluster-period of each outcome, respectively.
LCRT_NestExch = function(design_matrix,treatment_effect,time_effect_matrix,icc,cac,n,C){
  theta= treatment_effect
  prevs = design_matrix*theta+time_effect_matrix
  rhoct = icc
  SS = nrow(design_matrix)
  TT = ncol(design_matrix)
  if(cac<=1){
    rhoc = icc*cac
    # sample:
    data = matrix(-1,nrow=SS*C,ncol=TT*n) # matrix with one row per cluster. Observations in each row are grouped into time periods
    for(seq_num in 1:SS){
      data[(seq_num-1)*C+(1:C),] = NestBin::rNestBin(means=prevs[seq_num,],rhoC=rhoc,rhoCT=rhoct,n=n,C=C,sample=T)
    }
  }else{
    stop("cac must be less than or equal to 1")
  }
  # work out labels
  C_labels = matrix(rep(1:(SS*C),each=TT*n),nrow=SS*C,ncol=TT*n,byrow=TRUE)
  CT_labels = matrix(rep(1:(SS*C*TT),each=n),nrow=SS*C,ncol=TT*n,byrow=TRUE) # CT labels are 1:C*TT, not some ordered pair like (cluster,time)
  T_labels = matrix(rep(rep(1:TT,each=n),times=(SS*C)),nrow=(SS*C),ncol=TT*n,byrow=TRUE)
  df = data.frame(outcome=matrix(t(data),nrow=SS*C*TT*n),
                  C_labels =matrix(t(C_labels), nrow=SS*C*TT*n),
                  T_labels =matrix(t(T_labels), nrow=SS*C*TT*n),
                  CT_labels=matrix(t(CT_labels),nrow=SS*C*TT*n))
  return(df)
}