#' Toy example for rNestBin()
library("NestBin")
set.seed(1)
prevalences = c(0.1,0.2,0.3) # 3 prevalences, one per subcluster
rhoC = 0.1 # cluster-level correlation
rhoCT = 0.2 # subcluster-level correlation
n = 2 # number of observations per subcluster
C = 5 # number of clusters
Y = rNestBin(means=prevalences,rhoC=rhoC,rhoCT=rhoCT,n=n,C=C)
# > Y
#       [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]    0    0    0    1    0    1
# [2,]    0    1    0    0    0    0
# [3,]    0    0    0    0    1    1
# [4,]    0    0    0    0    0    0
# [5,]    0    0    0    0    0    0

p1 = 0.15
p2 = 0.126
prevalences2 = c(p1,p2)
rhoC = 0.025 # cluster-level correlation
rhoCT = 0.035 # subcluster-level correlation
n = 310
C = 1000
num_reps = 1000
mean_mat = matrix(NA,num_reps,n*length(prevalences2))
corr_vec_CT = rep(NA,num_reps)
corr_vec_C = rep(NA,num_reps)
for (i in 1:num_reps){
  Y = rNestBin(means=prevalences2,rhoC=rhoC,rhoCT=rhoCT,n=n,C=C)
  mean_mat[i,] = colSums(Y)/C
  corr_vec_CT[i] = cor(Y[,1],Y[,2])
  corr_vec_C[i] = cor(Y[,1],Y[,311])
}
# plotting results
old_par = par(no.readonly = TRUE)
par(mar=c(5,4,1,2)) # reduce empty space where there is no title
hist(mean_mat[,1],main=NA,xlab="Mean of first observation",mar=c(5, 4, 2, 2) + 0.1)
abline(v=p1,col='blue')
abline(v=mean(mean_mat[,1]),col='black')
hist(mean_mat[,311],main="",xlab="Mean of 311th observation")
abline(v=p2,col='blue')
abline(v=mean(mean_mat[,311]),col='black')
hist(corr_vec_C,main="",xlab="Cluster-level correlation, rhoC")
abline(v=rhoC,col='blue')
abline(v=mean(corr_vec_C),col='black')
hist(corr_vec_CT,main="",xlab="Cluster-period-level correlation, rhoCT")
abline(v=rhoCT,col='blue')
abline(v=mean(corr_vec_CT),col='black')
par(old_par)