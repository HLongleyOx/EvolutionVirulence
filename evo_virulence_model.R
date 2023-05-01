
#Install libraries
library(tidyverse)
library(deSolve)
library(abind)
library(pracma)

#Set seed
set.seed(1004)

#Load functions
source("tools/analytical_sol_heteropop.R")
source("tools/between_host_model.R")
source("tools/bh_equilibria.R")
source("tools/force_of_infection_calculation.R")
source("tools/hill_functions.R")
source("tools/many_mutation_model_calc.R")
source("tools/quasispecies_model_equilibrium_full.R")
source("tools/quasispecies_model_full.R")
source("tools/viralload_over_time.R")

#set fixed parameter values
max_dur_years <- (0.25 + hill_down_function(10^2) + 0.75) #max duration in yrs 
max_dur <- max_dur_years*365
num_strains_A=150
s_A=10^-4
num_strains_B=35
s_B=4*10^-3
num_strains_C=10
s_C=10^-1
mu=4*10^-5


vl_A = 7-(0:(num_strains_A))*0.033 #spvl approximatley goes between 
vl_B = 7-(0:(num_strains_B))*0.14
vl_C = 7-(0:(num_strains_C))*0.5

fullmodelA <- quasispecies_model_full(m=num_strains_A,M=num_strains_A+1,s=s_A, max_dur=20000, dt=1, mu=mu)
fullmodelB <- quasispecies_model_full(m=num_strains_B,M=num_strains_B+1,s=s_B, max_dur=max_dur, dt=1, mu=mu)
fullmodelC <- quasispecies_model_full(m=num_strains_C,M=num_strains_C+1,s=s_C, max_dur=max_dur, dt=1, mu=mu)



equilA <- cbind(quasispecies_model_equilibrium_full(m=num_strains_A,M=num_strains_A+1,s=s_A, mu=mu), spvl=vl_A)
equilB <- cbind(quasispecies_model_equilibrium_full(m=num_strains_B,M=num_strains_B+1,s=s_B, mu=mu), spvl=vl_B)
equilC <- cbind(quasispecies_model_equilibrium_full(m=num_strains_C,M=num_strains_C+1,s=s_C, mu=mu), spvl=vl_C)
equil <- rbind(equilA, equilB, equilC) #dist at equilibrium

vl_time_A <- viralload_over_time(num_strains_A, fullmodelA, vl_A)
vl_time_B <- viralload_over_time(num_strains_B, fullmodelB, vl_B)
vl_time_C <- viralload_over_time(num_strains_C, fullmodelC, vl_C)

#extended version of the model, 30,000 days: 
FOI_A = force_of_infection_calculation(time=max_dur, spvl=vl_A, muts_freq=fullmodelA, E=0, num_strains_A+1 )
FOI_B = force_of_infection_calculation(time=max_dur, spvl=vl_B, muts_freq=fullmodelB, E=0, num_strains_B+1)
FOI_C = force_of_infection_calculation(time=max_dur, spvl=vl_C, muts_freq=fullmodelC, E=0, num_strains_C+1)

durations_A <- FOI_A[[2]]
durations_B <- FOI_B[[2]]
durations_C <- FOI_C[[2]]

spvl_overTimeA = FOI_A[[3]]
spvl_overTimeB = FOI_B[[3]]
spvl_overTimeC = FOI_C[[3]]

spvl_A = c()
spvl_B = c()
spvl_C = c()

for (i in 1:(num_strains_A+1)){
  d = durations_A[i]*365
  spvl_A[i] = mean(spvl_overTimeA[1:d,i])
}

for (i in 1:(num_strains_B+1)){
  d = durations_B[i]*365
  spvl_B[i] = mean(spvl_overTimeB[1:d,i])
}

for (i in 1:(num_strains_C+1)){
  d = durations_C[i]*365
  spvl_C[i] = mean(spvl_overTimeC[1:d,i])
}


#between host model at equilibria: 
equil_prev_A <- bh_equilbria(FOI_A, durations_A,  spvl_A,num_strains_A+1, "A")
equil_prev_B <- bh_equilbria(FOI_B, durations_B, spvl_B,num_strains_B+1, "B")
equil_prev_C <- bh_equilbria(FOI_C, durations_C, spvl_C,num_strains_C+1, "C") #R0 too small
equil_prev <- rbind(equil_prev_A[[1]], equil_prev_B[[1]], equil_prev_C[[1]])


#between host model 
time_taken = 50 #time step
dt=time_taken/365
AICint = max_dur_years  #Duration of initial conditions: the max duration of the infection
iT0 = ceil(AICint/dt) #Time in days
num_host_types = 1
inds <- c(TRUE,rep(FALSE, time_taken-1))
fA = FOI_A[[1]][,,inds]
fB = FOI_B[[1]][,,inds]
fC = FOI_C[[1]][,,inds]
modA = list()
modB = list()
modC = list()

# #start with high viral load/average viral load/low viral load
initialVL_A = data.frame(sample(1:50, 10), initialVL_A = sample(50:90, 10), initialVL_A = sample(110:150, 10))
initialVL_B= data.frame(sample(1:10, 10), initialVL_A = sample(11:23, 10), initialVL_A = sample(26:36, 10))
initialVL_C = data.frame(sample(1:3, 3), initialVL_A = sample(4:7, 3), initialVL_A = sample(8:10, 3))

for (i in 1:3){
 modA[[i]] = between_host_model(f=fA,
                                 tmax=150,
                                 tstep=time_taken,
                                 max_dur=max_dur,
                                 initCons=initialVL_A[,i],
                                 num_strains=num_strains_A+1,
                                 num_hosts=1,
                                 durations_i=durations_A,
                                 host_dist=1,
                                 transmissibility=1)
}

for (i in 1:3){
  modB[[i]] = between_host_model(f=fB,
                                 tmax=100,
                                 tstep=time_taken,
                                 max_dur=max_dur,
                                 initCons=initialVL_B[,i],
                                 num_strains=num_strains_B+1,
                                 num_hosts=1,
                                 durations_i=durations_B,
                                 host_dist=1,
                                 transmissibility=1)
}

for (i in 1:3){
  modC[[i]] = between_host_model(f=fC,
                                 tmax=100,
                                 tstep=time_taken,
                                 max_dur=max_dur,
                                 initCons=initialVL_C[,i],
                                 num_strains=num_strains_C+1,
                                 num_hosts=1,
                                 durations_i=durations_C,
                                 host_dist=1,
                                 transmissibility=1)
}


host_effect <- analytical_sol_heteropop(50, num_strains_A, scenarioExtra_fullmodelA,spvl_A,E=0.5)


