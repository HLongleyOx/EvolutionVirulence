
#analytical solution to the between host model when we assume host hetrogeneity 
#host_type = number of host types 
#num strains = number of infection types
#model = within host model output
#spvl = spvl of virus type
#E = host effect size
#returns dataframe of of total prevalence, prevalence by mutation type, spvl distribution 

analytical_sol_heteropop <- function(host_types, num_strains, model,spvl,E) {
  
  host_effect = seq(-E, E, length.out = host_types)
  FOI_HH = lapply(
    host_effect,
    force_of_infection_calculation,
    time = max_dur,
    spvl = spvl,
    muts_freq = model,
    num_strains = num_strains + 1
  )
  
  durations_HH = array(dim = c(num_strains + 1, host_types))
  spvl_HH = array(dim = c(max_dur, num_strains + 1, host_types))
  Km_HH = array(dim = c(num_strains + 1, num_strains + 1, host_types))
  
  for (i in 1:host_types) {
    durations_HH[, i] = FOI_HH[[i]][[2]]
    spvl_HH[, , i] = FOI_HH[[i]][[3]]
    Km_HH[, , i] = FOI_HH[[i]][[4]]
  }

  
  KM_HH = list()
  durations_hetro = array(dim = c(num_strains + 1, host_types))
  spvl_hetro = array(dim = c(num_strains + 1, host_types))
  
  for (n in 1:host_types) {
    #find average viral load during infection
    KM_HH[[n]] = Km_HH[, , n]
    durations_hetro[, n] = durations_HH[, n]
    spvl_HH_overtime = spvl_HH[, , n]
    for (i in 1:(num_strains + 1)) {
      du = durations_hetro[i, n] * 365
      spvl_hetro[i, n] = mean(spvl_HH_overtime[1:du, i])
    }
  }
  
  spvl_hetro <-
    array(spvl_hetro, dim = c(host_types * (num_strains + 1), 1))
  uni = 1 / host_types #Assume uniform for now
  
  host_dist = rep(1 / host_types, host_types)
  for (h in 1:host_types) {
    KM = KM_HH[[h]]
    for (k in 1:host_types) {
      if (k == 1) {
        KM_h = uni * KM
      }
      else{
        KM_h = rbind(KM_h, host_dist[h] * KM)
      }
    }
    if (h == 1) {
      KM_combined = KM_h
    }
    else{
      KM_combined = cbind(KM_combined, KM_h)
    }
  }
  
  durations_H = array(durations_hetro, dim = c(host_types * (num_strains +
                                                               1), 1))
  d = 0.02 #natural mortality rate
  B = 200 #rate at which individuals enter susceptible population
  natural_mortality_last_HH <- exp(-(d) * durations_H)
  eig_km <- eigen(KM_combined)
  I = Re((eig_km$vec[, 1]))
  R_0 <- max(Re(eig_km$values))
  I = I / sum(I)
  I_strat_norm = I
  I_strat_adjust <-
    sum(I_strat_norm * natural_mortality_last_HH) #Sum of incidence adjusted for natural mortality
  I_tot <- B * (R_0 - 1) / (R_0 - I_strat_adjust)
  I_strat <- I_tot * I_strat_norm
  N_equil <- (B - I_tot * I_strat_adjust) / d #Size of population
  S_equil <- N_equil / R_0  #Number susceptible
  s_equil <- 1 / R_0 #fraction susceptible
  Ptot_equil <- N_equil - S_equil  #Prevalence
  ptot_equil <- 1 - s_equil
  Pbar_equil <-
    I_strat_norm * durations_H / sum(I_strat_norm * durations_H) #strst. normalised prev.
  P_equil <- Ptot_equil * Pbar_equil
  
  host_effect_bh <- data.frame(spvl=spvl_hetro, prevalence=P_equil,N=Ptot_equil)
  
  return(host_effect_bh)
}