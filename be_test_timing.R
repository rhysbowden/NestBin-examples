# be_test_timing: tries to compare time taken for 2 methods
# v4: cleaned up for use as supplementary script in paper submission
library('tictoc') # for timing
library('bindata') # one version of the Emrich and Piedmonte dichotomisation method
library('mvtBinaryEP') # another version of the Emrich and Piedmonte dichotomisation method
library('tidyverse') # for graphing the results
set.seed(4)
# structure settings
num_reps=1 # number of times to call each method for each combination of parameters
match_total = 1 # if 1, design 'none' has twice as many clusters per sequence to make up for having half as many sequences.
plot_results=FALSE
save_plot=FALSE
inflation_factor = 1000 # run faster methods inflation_factor times more, in order to get a better average time measurement
                        # impractical to run slower methods this many times

var_names = c('mu','rhoct','cac','cc','TT','nn','design')
# --- model/design parameters. Each one can be a vector, the simulations will run over all combinations of these parameters ---
mu_vals = c(0.4)
rhoct_vals = c(0.05)
cac_vals = c(0.5)
cc_vals = c(1)
nn_vals = c(20,40,60,80,100)
design_vals = c('cxo')
TT_vals = 2# 8 # must be even
stopifnot(all(TT_vals%%2==0))
param_grid = expand.grid(mu_vals,rhoct_vals,cac_vals,cc_vals,TT_vals,nn_vals,design_vals)
names(param_grid) = var_names
num_param_tuples = dim(param_grid)[1]

# observation parameters
theta = 0.05
num_methods=4
param_times = matrix(NA,nrow=num_param_tuples,ncol=num_methods)
method_times = rep(NA,num_methods)
for(method_ind in c(1,2)){
  tic(0)
  run_nestbin=0
  run_rmvbin=0
  run_rmvbin2=0
  run_mbe = 0
  if(method_ind == 1){
    run_rmvbin=1
  }
  if(method_ind == 2){
    run_nestbin=1
  }
  if(method_ind ==3){
    run_rmvbin2=1  # don't prestore results for rmvbin, do each rep separately
  }
  if(method_ind==4){
    run_mbe = 1 # runs mvtBinaryEP package version of Emrich and Piedmonte
  }
  for(param_ind in 1:num_param_tuples){
    cat(param_ind)
    mu = param_grid[param_ind,'mu']
    rhoct = param_grid[param_ind,'rhoct']
    cac = param_grid[param_ind,'cac']
    rhoc = rhoct*cac
    cc = param_grid[param_ind,'cc']
    nn = param_grid[param_ind,'nn']
    TT = param_grid[param_ind,'TT']
    design_name = param_grid[param_ind,'design']
    
    # design parameters
    beta = rep(c(0.05,0),TT/2)#c(0.1,0,0.1,0,0.1) # time fixed effects, must be TT long
    if(design_name =='none'){
      SS = 1 # number of sequences
      sequences = matrix(0,nrow=1,ncol=nn*TT) # treatment sequences, one per row
      sequence_reps = rep((1+match_total)*cc,SS) # number of clusters in each sequence
    }
    if(design_name =='cxo'){
      stopifnot(TT%%2 == 0)
      SS = 2 # number of sequences
      sequences = rbind(rep(c(0,1),each=nn*TT/2),rep(c(1,0),each=nn*TT/2))# treatment sequences, one per row
      sequence_reps = rep(cc,SS) # number of clusters in each sequence
    }
    if(design_name =='pa'){
      SS = 2 # number of sequences
      sequences = rbind(rep(0,nn*TT),rep(1,nn*TT)) # treatment sequences, one per row
      sequence_reps = rep(cc,SS) # number of clusters in each sequence
    }
    CC = sum(sequence_reps) # total number of clusters
    Ncells = CC*TT
    NN = Ncells*nn
    
    # vector versions of fixed effects, for one cluster
    beta_vecC = rep(beta,each=nn)
    mu_vecC = rep(mu,TT*nn)
    marg_probs = matrix(-1,nrow =SS,ncol=nn*TT) # marginal probabilities for each observation in a cluster, one row per treatment sequence
    for(ss in 1:SS){
      marg_probs[ss,] = mu_vecC + beta_vecC + theta*sequences[ss,]
    }
    
    # nested cluster correlation structure
    time_labelsC = rep(1:TT,each=nn)
    time_labelsCMhorz = matrix(rep(time_labelsC,each=TT*nn),nrow=TT*nn,ncol=TT*nn)
    time_labelsCMvert = matrix(rep(time_labelsC,TT*nn),nrow=TT*nn,ncol=TT*nn)
    corr_matCT = matrix(rhoc,nrow=nn*TT,ncol=nn*TT)
    corr_matCT[time_labelsCMhorz==time_labelsCMvert] = rhoct
    diag(corr_matCT)=1
    
    # run methods once to get over initialization problems
    if(run_nestbin){
      ccs = sequence_reps[1]
      PPs = marg_probs[1,]
      temp =rNestBin(means=PPs[(1:TT)*nn], rhoC=rhoc, rhoCT=rhoct, n = nn, C = ccs, sample = TRUE)
    }
    if(run_rmvbin|run_rmvbin2){
      temp = rmvbin(n=sequence_reps[ss],margprob=marg_probs[ss,],bincorr=corr_matCT)
    }
    if(run_mbe){
      temp = ep(mu=marg_probs[ss,],R=corr_matCT)
    }
    
    # start time!
    tic(param_ind)
    if(run_rmvbin){
      YY2m_total = list()
      for(ss in 1:SS){
        YY2m_total[[ss]] = rmvbin(n=sequence_reps[ss]*num_reps,margprob=marg_probs[ss,],bincorr=corr_matCT)
      }
    }
    
    # main loop
    for(repnum in 1:(num_reps*inflation_factor^(method_ind==2|method_ind==5))){
      # sampling
      # means
      YY2m_part = list()
      for(ss in 1:SS){
        ccs = sequence_reps[ss]
        PPs = marg_probs[ss,]
        if(run_rmvbin){
          YY2m_part[[ss]] = YY2m_total[[ss]][(repnum-1)*ccs+(1:ccs),]
        }
        if(run_rmvbin2){
          YY2m_part[[ss]] = rmvbin(n=sequence_reps[ss],margprob=marg_probs[ss,],bincorr=corr_matCT)
        }
        if(run_nestbin){
          YY2m_part[[ss]] = rNestBin(means=PPs[(1:TT)*nn], rhoC=rhoc, rhoCT=rhoct, n = nn, C = ccs, sample = TRUE)
        }
        if(run_mbe){
          YY2m_part[[ss]] = ep(mu=marg_probs[ss,],R=corr_matCT) # could store isd for each sequence to speed up future reps
        }
      }
      YY2m = do.call(rbind,YY2m_part)
      YY2 = as.vector(YY2m)
    }
    temp_time = toc(param_ind)
    param_times[param_ind,method_ind] = unname(temp_time$toc-temp_time$tic)
  }
  temp_time2 = toc(0)
  method_times[method_ind] = unname(temp_time2$toc-temp_time2$tic)
}

# plotting:
if(plot_results){
  pt= data.frame(param_times)
  pt[,2] = pt[,2]/inflation_factor
  pt[,3] = param_grid$nn
  names(pt) = c('Emrich and Piedmonte','NestBin','n')
  ptL = pivot_longer(data=pt,cols=c('Emrich and Piedmonte','NestBin'),names_to='method',values_to='time')
  plot2 = ggplot(ptL)+
    geom_line(mapping=aes(x=n,y=time,colour=method,linetype=method))+
    scale_y_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000),
                  labels=c('0.00001','0.0001','0.001','0.01','0.1','1','10','100','1000','10000') )+
    labs(x='observations per subcluster',y='time per call (seconds)')+
    theme_bw()+
    theme(legend.position = "none")
  print(plot2)
  if(save_plot){
    ggsave(paste("../timing_plot",".pdf",sep=""), plot=plot2, height=10, width=20, units="cm", dpi=600)
  }
}