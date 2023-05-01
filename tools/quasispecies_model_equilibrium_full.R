# Equilibrium solution of the quasispecies equation
# Return dataframe with frequency of each virus type at the equilibrium solution
# m = max number of mutations 
# M = number of virus types (m+1)
# s = ftiness cost
# mu = mutation rate

quasispecies_model_equilibrium_full <-  function(m, M, s, mu = mu) {
  mu = mu

    n = m + 1
  ivector <- 0:m
  s <- s
  A <- (1 - s) ^ ivector

  Q = matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    Q[i, i] <- (1 - mu) ^ (m)  #Error-free replication
    Q[i, i - 1] <-
      mu * (m - (i - 2)) * (1 - mu) ^ (m - 1) #Acquiring one further mutation
    if (i < n) {
      Q[i, i + 1] <-
        mu * (i) * (1 - mu) ^ (m - 1) #Back mutation at one of the already mutated sites
    }
  }
  
  W = Q %*% diag(A)
  
  W = W[1:M, 1:M]
  A = A[1:M]
  
  eigenW <- eigen(W)
  vec <- eigenW$vectors[, 1]
  vec <- abs(Re(vec))
  vec2 <- vec / sum(vec)
  mutations = 0:(M - 1)
  mean <- sum(vec2 * mutations)
  sd <- sum((mutations - mean) ^ 2 * vec2)
  df <-
    data.frame(
      "equilibrium" = vec2,
      mutations,
      mean = mean,
      sd = sd,
      h = vec
    )
  return(data.frame(df))
} 

