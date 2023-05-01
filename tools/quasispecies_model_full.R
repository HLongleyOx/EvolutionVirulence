
# Set up quasispecies model for every possible starting conditions 
# Outputs 3 dimensional matrix with solution of quasispecies equation for every
# infection type
# mu = mutation rate
# m = max number of mutations
# s = fitness coefficient
# M = num strains 
# max dur = max duration (generations)
# dt = time step

quasispecies_model_full <- function(m,M,s, max_dur, dt, mu){
  
  mu=mu
  n=m+1
  ivector <- 0:m
  s <- s
  A <- (1-s)^ivector
  Q = matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n){
    Q[i,i] <- (1-mu)^(m)  #Error-free replication
    Q[i,i-1] <- mu*(m-(i-2))*(1-mu)^(m-1) #Acquiring one further mutation
    if (i<n) {
      Q[i,i+1] <- mu*(i)*(1-mu)^(m-1) #Back mutation at one of the already mutated sites
    }
  }
  W = Q %*% diag(A)
  
  Yi_mat <- diag(rep(1,M))
  many_mutation_model<-apply(Yi_mat, MARGIN=1, FUN=many_mutation_model_calc,W=W,mu=mu,
                             m=m, M=M,dt=dt, max_dur=max_dur)
  
  finalmodel <- list(many_mutation_model)
  finalmodel <- finalmodel[[1]]
  virus_freq=array(dim=c(M,M,max_dur))
  
  for (k in 1:max_dur){
    for (j in 1:M){
      df <- finalmodel[[j]][k,]   #frequency of type j infection at time i
      virus_freq[,j,k] <- as.numeric(df[1,-1])
    }
  }
  return(virus_freq)
} 
