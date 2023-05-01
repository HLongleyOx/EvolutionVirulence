
#Find viral load over time 
#num_strains = number of  virus types
#mod output describes the frequency of virus type j at time t of a type i  infection
#spvl is the viral load of each virus type

viralload_over_time <- function(num_strains, mod_output, spvl){
  spvl_time <- data.frame(spvl=c(), time=c(), start=c())
  for (i in 1:(num_strains+1)){
    df <- mod_output[,i,]
    df_spvl <- apply(df,2,function(x,spvl){sum(x*spvl)}, spvl=spvl)
    spvl_time_temp <- data.frame(spvl=df_spvl, time=seq(1, length(df_spvl)), start=i)
    spvl_time <- rbind(spvl_time, spvl_time_temp)
  }
  return(spvl_time)
}
