# be_test_search.R
# searches for a prevalence that allows each method to work for a given rhoc,rhoct
# the methods are: 
# rNestBin from the NestBin package (Bowden2022)
# CBern from the CorBin package (Jiang2021)
# the conditional linear family method from Qaqish (with various choices for the ordering of the binary RVs)

library("CorBin")
library('tidyverse')

# qaqish_feasR(rhoC,rhoCT,beta,mm,num_samples,tol=0,include_id=TRUE)
#
# Finds whether a given mean and nested exchangeable correlation structure for a
# sequence of binary random variables can be represented by B. F. Qaqish's 
# Conditional Linear Family (CLF) of models.
#
# The method for testing this is given in 
# Qaqish (2003), A family of multivariate binary distributions for simulating
# correlated binary variables with specified marginal means
# and correlations.
# Whether or not a given structure can be represented depends on the ordering of
# the binary random variables to be simulated, and so this function tries 
# several permutations of those variables and declares the structure to be 
# feasible if any one of those permutations allows for representation using
# a CLF model. In this case, TRUE is returned, otherwise FALSE is returned.
#
# qaqish_feasR(rhoC,rhoCT,beta,mm,num_samples,tol,include_id=TRUE)
#
# inputs:
# rhoC: cluster-level correlation
# rhoCT: subcluster-level correlation, rhoCT>rhoC
# beta: the vector of prevalences, one for each subcluster. len(beta) determines
#       the number of subclusters per cluster.
# mm: the number of observations in each subcluster.
# num_samples: the number of random permutations to try.
# tol: a tolerance for when testing that critical values fall between 0 and 1.
# include_id: whether or not the identity permutation is guaranteed to be 
#             included.
qaqish_feasR = function(rhoC,rhoCT,beta,mm,num_samples,tol=0,include_id=TRUE){
  # Example parameters:
  #beta = c(0.1,0.2,0.4) # beta(i) is prevalence in period i
  #mm = 10 # number per cluster period
  #rhoCT = 0.2
  #rhoC = 0.1
  is_feas = rep(FALSE,num_samples)
  for(sample_num in 1:num_samples){
    TT=length(beta)
    NN = mm*TT
    beta_i = rep(beta,each=mm)
    cor_i = rhoC*matrix(1,mm*TT,mm*TT)+(rhoCT-rhoC)*diag(1,TT)%x%matrix(1,mm,mm)+(1-rhoCT)*diag(1,mm*TT)
    sigma_vec = sqrt(beta_i*(1-beta_i))
    cov_i = diag(sigma_vec)%*%cor_i%*%diag(sigma_vec)
    
    # random permutation of observations:
    if(i==1){
      permut = 1:ncol(cov_i)
    }else{
      permut = sample(ncol(cov_i))
    }
    Rcov_i = cov_i[permut,permut]
    Rbeta_i = beta_i[permut]
    
    bij = list()
    for(i in 2:NN){
      bij[[i]] = solve(Rcov_i[1:(i-1),1:(i-1)])%*%Rcov_i[1:(i-1),i]
    }
    max_values = rep(0,NN)
    min_values = rep(0,NN)
    for(i in 2:NN){
      max_values[i] = Rbeta_i[i]+sum(bij[[i]]*(1-Rbeta_i[1:(i-1)]))
      min_values[i] = Rbeta_i[i]-sum(bij[[i]]*Rbeta_i[1:(i-1)])
    }
    is_feas[sample_num] = all(max_values<=(1+tol))&all(min_values>=(0-tol))
  }
  return(any(is_feas))
}

# determines whether a given method can simulate with a certain set of parameters
method_fails = function(rhoct,mu,beta_max,method,nn,TT){
  match_total=0
  run_bef = F
  run_cBern = F
  run_qaq = F
  run_qaqR = F
  mfail=0
  if(method=='B'){
    run_bef = T
  }
  if(method=='C'){
    run_cBern = T
  }
  if(method=='Q'){
    run_qaq = T 
  }
  if(method=='rQ'){
    run_qaqR = T
  }
  #mu = 0.2
  #rhoct = 0.024
  cac = 0.8
  rhoc = rhoct*cac
  cc = 1
  theta = 0.2
  design_name = 'none'
  
  # design parameters 
  beta = rep(c(beta_max,0),TT/2)#c(0.1,0,0.1,0,0.1) # time fixed effects, must be TT long
  if(TT%%2==1){
    beta = c(beta,beta_max)
  }
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
  
  YY2m_partC = list()
  YY2m_partB = list()
  for(ss in 1:SS){
    ccs = sequence_reps[ss]
    PPs = marg_probs[ss,]
    if(run_cBern){
      # CorBin, cBern version
      # generate list of diagonals of correlation matrix:
      diag_list = list()
      mm = dim(corr_matCT)[1]
      for(h in 1:(mm-1)){
        diag_list[[h]] = c(corr_matCT[col(corr_matCT)==(row(corr_matCT)+h)])
      }
      # sample
      YY2m_partC[[ss]] = cBern(n=ccs,p=PPs,rho=diag_list,type="General",k=mm-1)
      if(!is.numeric(YY2m_partC[[ss]])){
        mfail=1
      }
    }
    if(run_bef){
      this_mfail = tryCatch({
        YY2m_partB[[ss]] = NestBin(means=PPs[(1:TT)*nn],rhoC=rhoc,rhoCT=rhoct,n=nn,C=ccs,sample=F)
        return(0)
      },
      warning=function(war){
        #converge_fail2[i]=1
        print(paste("MY_WARNING:  ",war))
        return(0)
      },
      error=function(err){
        return(1)
      })
      mfail = mfail+this_mfail
    }
    if(run_qaq){
      mfail = mfail+(1-qaqish_feas(rhoC=rhoc,rhoCT=rhoct,beta=PPs[(1:TT)*nn],mm=nn,tol=1E-7))
    }
    if(run_qaqR){
      mfail = mfail+(1-qaqish_feasR(rhoC=rhoc,rhoCT=rhoct,beta=PPs[(1:TT)*nn],mm=nn,num_samples=20,tol=1E-7))
    }
  }
  return(mfail)
} ############### function end ###############

# find feasible region
mu_val=0.1
beta_max=0.2
n_vals = c(10)#,20,30,40,50)
T_vals = c(5)#,10,15,20)
max_feasC = matrix(0,nrow=length(n_vals),ncol=length(T_vals))
max_feasB = matrix(0,nrow=length(n_vals),ncol=length(T_vals))
max_feasQ = matrix(0,nrow=length(n_vals),ncol=length(T_vals))
max_feasQB = matrix(0,nrow=length(n_vals),ncol=length(T_vals))
max_feasRQ = matrix(0,nrow=length(n_vals),ncol=length(T_vals))
for(i in 1:length(n_vals)){
  for(j in 1:length(T_vals)){
    nn = n_vals[i]
    TT = T_vals[j]
    testQ = function(rhoct){method_fails(rhoct=rhoct,mu=mu_val,beta_max=beta_max,method='Q',nn=nn,TT=TT)-0.5}
    testRQ = function(rhoct){method_fails(rhoct=rhoct,mu=mu_val,beta_max=beta_max,method='rQ',nn=nn,TT=TT)-0.5}
    
    if(testQ(0.0001)>0){
      max_feasQ[i,j]=0
    }
    else if(testQ(0.9999)==0){
      max_feasQ[i,j]=0.9999
    }else{
      max_feasQ[i,j] = uniroot(testQ,c(0.0001,0.9999))$root
    }

    if(testRQ(0.0001)>0){
      max_feasRQ[i,j]=0
    }
    else if(testRQ(0.9999)==0){
      max_feasRQ[i,j]=0.9999
    }else{
      max_feasRQ[i,j] = uniroot(testRQ,c(0.0001,0.9999))$root
    }
  }
}

# facet plot:
#long_max_feasB = matrix(max_feasB,ncol=1)
#long_max_feasC = matrix(max_feasC,ncol=1)
long_max_feasQ = matrix(max_feasQ,ncol=1)
long_max_feasRQ = matrix(max_feasRQ,ncol=1)
#long_beta_max_vals = rep(beta_max_vals,each=nrow(max_feasB))
#long_mu_vals = rep(mu_vals,ncol(max_feasB))
long_nn_vals = rep(n_vals,times=length(T_vals))
long_TT_vals = rep(T_vals,each=length(n_vals))
#df1 = data.frame(min_val= long_mu_vals,range=long_beta_max_vals,BEBin = long_max_feasB,CorBin=long_max_feasC,Qaqish=long_max_feasQ,RQaqish=long_max_feasRQ)
#df1 = data.frame(nn_val= long_nn_vals,T=long_TT_vals,BEBin = long_max_feasB,Qaqish=long_max_feasQ,Jiang=long_max_feasC)
df1 = data.frame(nn_val= long_nn_vals,T=long_TT_vals,Qaqish=long_max_feasQ,RQaqish=long_max_feasRQ)
df2 = pivot_longer(df1,cols=c('Qaqish','RQaqish'),names_to = 'method',values_to='max_rho')
#df2$method = factor(df2$method,levels=c("BEBin","CorBin","CLF"))
#df2=df2[df2$method!='BEBin',]
df2$method = factor(df2$method,levels=c("Qaqish","RQaqish"))
plot1 = ggplot(df2,mapping=aes(x=nn_val,y=max_rho,group=method))+
  geom_line(aes(color=method,linetype=method))+
  facet_wrap(~T,labeller=label_both)+
  ylab('Maximum correlation')+
  xlab('N')+
  theme_bw()+
  theme(legend.position = "none")
#plot1+labs(x='Maximum correlation')
print(plot1)
if(0){
ggsave(paste("../rplots/correlation_rangesBQC_nn10-50_TT5-10_no_legend.pdf"), plot=plot1, height=10, width=21, units="cm", dpi=600)
}