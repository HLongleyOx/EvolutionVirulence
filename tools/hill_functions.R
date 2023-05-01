# hill functions to determine infectivity and duration
# Parameterised from Fraser et al 2007 

hill_up_function <- function(x,Cmax=0.317,C50=13938,k=1.02){
  output = (Cmax * x^k)/( x^k + C50^k)
  return(output)
}

#For duration
hill_down_function <- function(x,Dmax=25.4,D50=3058,Dk=0.41){
  output = (Dmax * D50^Dk)/(x^Dk + D50^Dk) 
  return(output)
}
