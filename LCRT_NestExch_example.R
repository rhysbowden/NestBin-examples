#' Toy example for LCRT_NestExch
#' 
#' Simulates multiple stepped wedge designs using LCRT_NestExch() and then plots
#'  the observed distribution of correlations between two observations in the 
#'  same cluster-period, and two observations in differing periods but the same 
#'  cluster.
source("LCRT_NestExch.R")
num_reps=1000
num_sequences = 2
num_periods = 3
corr_vec_C = rep(NA,num_reps)
corr_vec_CT = rep(NA,num_reps)
for(j in 1:num_reps){
  num_trials = 1000
  design_matrix = matrix(c(0,1,1,0,0,1),nrow=num_sequences,ncol=num_periods,byrow=T)
  treatment_effect = 0.05 # additive treatment effect
  time_effect_matrix = matrix(rep(c(0.11,0.15,0.17),2),nrow=num_sequences,ncol=num_periods,byrow=T) # additive time effects, the (i,j)th element of which is the time period for cell
  icc=0.1 # example within-cluster-period correlation
  cac=0.6 # cluster auto-correlation, ratio of within- to between-cluster-period correlations 
  n=2 # number of observations per cluster-period cell
  C=3 # number of clusters per treatment sequence
  outcomes = matrix(0,nrow=num_trials,ncol=nrow(design_matrix)*ncol(design_matrix)*C*n) # one row per trial
  for(i in 1:num_trials){
    df = LCRT_NestExch(design_matrix,treatment_effect,time_effect_matrix,icc,cac,n,C)
    outcomes[i,] = df$outcome
  }
  corr_vec_CT[j] = cor(outcomes[,25],outcomes[,26]) # store the correlation of two particular observations in the same cluster-period cell over all the num_trials trials
  corr_vec_C[j] = cor(outcomes[,22],outcomes[,24]) # store the correlation of two particular observations in the same cluster but different cluster-period cells
}

old_par = par()
par(mar=c(5,4,1,2)) # reduce empty space where there is no title

hist(corr_vec_C,main="",xlab="Cluster-level correlation, rhoC")
abline(v=icc*cac,col='blue')
abline(v=mean(corr_vec_C),col='black')

hist(corr_vec_CT,main="",xlab="Cluster-period-level correlation, rhoCT")
abline(v=icc,col='blue')
abline(v=mean(corr_vec_CT),col='black')

par(old_par)