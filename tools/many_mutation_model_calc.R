# Function to numerically solve ODEs given initial conditions. Specifically, 
# solves the quasispecies equations and outputs a dataframe with results.
# Returns dataframe of frequencies of each virus type every generation 
# W=reproduction-mutation matrix
# Yi=initial frequencies of each virus type
# mu = mutation rate
# m = max number of mutations
# M = num strains 
# max dur = max duration (generations)
# dt = time step


many_mutation_model_calc <-
  function(W, Yi, mu, m, M, max_dur, dt) {
    SI <- function(time, state, params) {
      with(as.list(c(time, state, params)), {
        Yi = state
        w_prod = sum(W %*% Yi)
        dYi = W %*% Yi - w_prod * Yi
        return(list(c(dYi)))
      })
    }
    inits = c(Yi)
    params = c(
      W = W,
      mu = mu,
      m = m,
      M = M
    )
    times = seq(0, max_dur, by = dt)
    finalmodel = ode(inits, times, SI, params, atol = 10 ^ (-15))
    final_model <- as.data.frame(finalmodel)
    return(final_model)
  }