#determine force of infection at every time point, based upon viral load 
#allows for host hetrogeneity 
#time = max number of generations
# spvl = vector of viral loads associated with each virus type
# muts_freq = matrix with frequency of each mutations type from time t=1 to time
# t = max_dur, for every possible infection type
# E = host effect 
# num_strains = number of mutation types
#returns a list of outputs:
#force_of_infection = the force of infection of every virus typ 
#durations_i = duration of infection type i
#spvl_i = spvl of infection type i at each time point 
#KM = matrix that describes expected number of onward infections of each type
#durations = duration of infection type i

force_of_infection_calculation <-
  
  function(time, spvl, muts_freq, E, num_strains) {
    spvl_i <- array(data = NA, dim = c(time, num_strains))
    infect_i <- array(data = NA, dim = c(time, num_strains))
    transmiss_i <- array(data = NA, dim = c(time, num_strains))
    
    for (i in 1:time) {
      spvl_i[i,] <- colSums(muts_freq[, , i] * spvl) + E
      infect_i[i,] <- sapply(10 ^ spvl_i[i,], hill_up_function)
      transmiss_i[i,] <- 1
    }
    
    acute_stage <- round(0.25 * 365)
    AIDS_stage <- round(0.75 * 365)
    infect_i[1:(acute_stage),] <- 0.76
    durations_i <- c()
    
    #remaining duration:
    for (i in 1:num_strains){
      durations_temp <- round((hill_down_function(10^spvl_i[1,i])+1)*365)
      if (durations_temp>time){
        durations_temp = max_dur
      }
      spvl_mean <- weighted.mean(spvl_i[(1:durations_temp),i], w=seq(durations_temp,1,-1)*(1/durations_temp))
      durations_i[i] <- hill_down_function(10^spvl_mean)+1
      if (durations_i[i]*365>time){
        durations_i[i] = (max_dur/365)
      }
    }
    
    durations_days <-  round(durations_i * 365)
    for (i in 1:num_strains) {
      dur <- durations_days[i]
      infect_i[(dur - AIDS_stage):time, i] = 0
    }
    
    #Create beta term for strain-specific infectivity profile
    beta <- array(NA, dim = c(num_strains, num_strains, time))
    for (k in 1:time) {
      beta[, , k] = muts_freq[, , k] * infect_i[k,] * 1
    }
    
    #Add time component
    d <- 0.02 / 365
    duration_profiles <-
      array(data = NA, dim = c(time, num_strains))
    
    for (i in 1:num_strains) {
      time_seq <- seq(0, (time - 1))
      duration_days <- durations_days[i]
      duration_profiles[, i] <-
        ifelse(time_seq > duration_days, 0,  exp(-d * time_seq))
    }
    
    #at time t since epidemic began, the force of infection of strain i in type j
    force_of_infection <-
      array(NA, dim = c(num_strains, num_strains, time))
    
    dt = 1 / 365
    durations <- round((durations_i) / dt)
    KM <- zeros(num_strains)
    
    for (j in 1:num_strains) {
      for (i in 1:num_strains) {
        asymp_dur = durations[j]
        if (asymp_dur == 0) {
          KM[i, j] = 0
          force_of_infection[i, j, ] = 0
        } else {
          force_of_infection[i, j, 1:asymp_dur] = ((infect_i[, j] * muts_freq[i, j, ] *
                                                      exp(-d * 0:(time -
                                                                    1))))[1:durations[j]]
          f1 =  force_of_infection[i, j, 1:asymp_dur]
          KM[i, j] <-  dt * trapz(y = f1, x = 1:(length(f1)))
        }
      }
    }
    return(list(force_of_infection, durations_i, spvl_i, KM, durations))
  }
