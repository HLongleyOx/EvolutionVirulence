#solution to between host equilibrium 
#FOI=force of infection matrix
#durations_i = durations of each infection
#spvl_i=spvl of infection type i
# m=number of virus types
#model = name of model
#returns the distribution of spvls at equilibrium, incidence and prevalence
#by infection type, and R0

bh_equilbria <- function(FOI, durations_i, spvl_i,m,model){
  
  d=0.02
  B=200
  dt = 1
  
  durations <- round((durations_i)/dt)
  duration_years <- durations_i
  natural_mortality_last <- exp(-(0.02*duration_years))
  
  Km <- FOI[[4]]
  eig_km <- eigen(Km)
  R_0 <- max(Re(eig_km$values))
  I_eig <- abs(Re(eig_km$vec[,1]))
  I_strat_norm <- I_eig/sum(I_eig)
  I_strat_adjust <- sum(I_strat_norm * natural_mortality_last) #Sum of incidence adjusted for natural mortality
  I_tot <- B *(R_0-1)/(R_0-I_strat_adjust)
  I_strat <- I_tot* I_strat_norm 
  N_equil <- (B-I_tot * I_strat_adjust)/d #Size of population
  S_equil <- N_equil/R_0  #Number susceptible
  s_equil <- 1/R_0 #fraction susceptible
  Ptot_equil <- N_equil-S_equil  #Prevalence
  ptot_equil <- 1-s_equil
  Pbar_equil <- I_strat_norm * duration_years / sum(I_strat_norm*durations_i) #strst. normalised prev.
  P_equil <- Ptot_equil * Pbar_equil
  spvl_equil <- sum(Pbar_equil * spvl_i)
  spvl_prevalence <- data.frame(spvl=spvl_i, prevalence = P_equil, model=model)
  spvl_prevalence$m=m
  return(list(spvl_prevalence, I_strat/sum(I_strat),P_equil, R_0))
}
