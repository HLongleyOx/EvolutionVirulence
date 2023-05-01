#Determines numerical solution to between host model
#Can account for different host types
#f = force of infection
#tmax = total duration for between-host epidemic
#tstep = time step
#max_dur=  max duration of infection 
#init cons = initial infection type distributon
#num strains = number of infection types
#num host = number of host types
#durations_i = duration of type i infection
#host_dist = distribution of host types in population
#all infection types are equally transmissible 
#returns list of dataframe with total N, S, I and P, 
#a matrix with stratified incidence at each time point
#a matrix with stratified prevalence at each time point 

between_host_model <-
  function(f,
           tmax,
           tstep,
           max_dur,
           initCons,
           num_strains,
           num_hosts,
           durations_i,
           host_dist,
           transmissibility=1 ) {
    
    tmin = 0
    tmax = tmax
    time_taken = tstep
    dt = time_taken / 365
    ICint = max_dur / 365  #Duration of initial conditions: the max duration of the infection
    AICint = max_dur / 365 #Time ins years
    iT0 = ceil(AICint / dt) #Time in days
    num_host_types = num_hosts
    
    t = seq(tmin - AICint, tmax, by = dt)
    tt = seq(0, ICint, by = dt)
    lt = length(t)
    
    B = 200  #birth rate
    d = 0.02 #deaths per cap per year
    natural_mortality <- exp(-(d) * tt)
    durations <- round((durations_i) / dt)
    duration_years <- durations_i
    natural_mortality_last <- exp(-(0.02 * (dt)) * durations)
    ltt = length(tt) #Max Time points in one infection
    
    I = array(0, dim = c(num_strains, lt)) #Incidence of each strain at every time point
    Ibar = array(0, dim = c(num_strains, lt))
    
    Itot = zeros(1, lt) #Total incidence across all strains
    
    P = array(0, dim = c(num_strains, lt)) #Incidence of each strain at every time point
    Pbar = array(0, dim = c(num_strains, lt))
    Ptot = zeros(1, lt) #Total incidence across all strains
    
    S = zeros(1, lt) #susceptibles at each time point
    N = zeros(1, lt) #Population size at each time point
    Y = array(0, dim = c(num_strains, lt))#Amount of virus in population, of each strain
    
    
    full_solution <- list()
    Incidence_stratified <- list()
    Prevalence_statified <- list()
    I0 <- zeros(num_strains, 1)
    
    I0[initCons] <- 1   #1 initial infection
    
    #Initial conditions
    N0 <- 10000 #10,000 initial susceptible
    S0 <- N0
    #Set ICs for first lot of time points
    for (iic in 1:iT0) {
      I[, iic] = I0
      Ibar[, iic] = I0 / sum(I0)
      Itot[iic] = sum(I0)
      N[iic] = N0
      S[iic] = S0
      P[, iic] = zeros(num_strains, 1)
      Pbar[, iic] = zeros(num_strains, 1)
      Ptot[iic] = 0
    }
    
    death_prob <- exp(-(d * (dt)) * tt)
    S[iT0] = N[iT0] - Itot[iT0] * dt * trapz(death_prob)
    
    FOI = array(0, dim = c(num_strains, lt))
    FOI_temp2 = zeros(num_strains)
    ltt = length(tt)
    
    #Loop through time points, then each strain, then each infecting strain
    for (it in (iT0 + 1):lt) {
      #each time point following the ICs
      discI = zeros(num_strains)
      for (i in 1:num_strains) {
        for (k in 1:num_host_types) {
          for (j in 1:num_strains) {
            FOI_temp = array(0, dim = c(num_strains, lt)) #Force of infection from type j infection in past
            FOI_temp[j, 1] = 0
            for (itt in 2:min(it, durations[j])) {
              FOI_temp[j, itt] =   f[i, j, itt] * I[j, (it - itt + 1)]  #force of infection * incidence(make 3d)
            }
            FOI_temp2[j, k] = dt * trapz(FOI_temp[j, 1:min(it, durations[j])]) #Number of new infections between time points
          }
        }
        for (k_rec in 1:num_host_types) {
          host_type_dist <- host_dist[k_rec]
          Y[i, it] = host_type_dist * sum(FOI_temp2) #Sum over infection & transmitted host types
          FOI[i, it] = transmissibility * Y[i, it]
          I[i, it] = S[it - 1] * FOI[i, it] / N[it - 1] #Incidence of type i infection in type k_rec host
          q = seq(it, (it - durations[i]), by = -1)
          P[i, it] = dt * trapz(I[i, q] *
                                  natural_mortality[1:(durations[i] + 1)])
          discI[i] = I[i, it - durations[i]] * natural_mortality_last[i]
        }
      }
      
      
      Itot[it] = sum(I[, it])   #Total incidence
      Ptot[it] = sum(P[, it])   #Total prevalence
      S[it] = N[it - 1] - Ptot[it]
      N[it] = N[it - 1] + dt * (B - (d) * N[it - 1] - sum(discI))  #Total populatio
      
      Ibar[, it] = I[, it] / Itot[it]
      Pbar[, it] = P[, it] / Ptot[it]
      
    }
    
    full_solution <-
      data.frame(
        Itot = t(Itot),
        Ptot = t(Ptot),
        S = t(S),
        N = t(N)
      )
    Incidence_stratified <- data.frame(Ibar)
    Prevalence_statified <- data.frame(Pbar)
    
    return(list(full_solution, Incidence_stratified, Prevalence_statified))
    
  }
