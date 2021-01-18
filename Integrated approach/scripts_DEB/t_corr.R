#'@export
t_corr <- function (p, Ti) {
  # this function calculates de correction factor for the input temperature. It
  # uses the vector t.pars yielded with the get_DEB_pars function. As input it
  # takes the list from get_DEB_pars and the temperature at which towards the par
  # or variable needs to be corrected (T)
  
  if (max(Ti<100)){ # then the units are likely not K
    Ti=Ti+273.15  # transform temperature to K
    warning('Ti was assumed to be in degrees C and converted to K')
  }
  
  if (length(p)>6){
  t.pars=p$t.pars  # extract vector with temperatures from p
  } else {t.pars=p}
  
  Ti = Ti
  
  if (length(t.pars)==2){  # only Arrhenius temperature, no boundaries are known
    T_ref = t.pars[1]
    T_A = t.pars[2]
    
    TC = exp(T_A/T_ref - T_A/Ti)
    
  }
  
  if (length(t.pars)==4){ # Arrhenius and the lower boundary temperartures with the changes in slope at such boundarie
    T_ref = t.pars[1]
    T_A = t.pars[2]
    T_L = t.pars[3]
    T_AL = t.pars[4]
      
    TC = exp(T_A/ T_ref - T_A/ Ti) *
      (1 + exp(T_AL/ T_ref - T_AL/ T_L))/
      (1 + exp(T_AL/ Ti - T_AL/ T_L))
  }
  
  if (length(t.pars)==6){ # Arrhenius and the boundary temperartures with the changes in slope at such boundaries
    T_ref = t.pars[1]
    T_A = t.pars[2]
    T_L = t.pars[3]
    T_H = t.pars[4]
    T_AL = t.pars[5]
    T_AH = t.pars[6]
    
    TC = exp(T_A/ T_ref - T_A/ Ti) *
      (1 + exp(T_AL/ T_ref - T_AL/ T_L) + exp(T_AH/ T_H - T_AH/ T_ref))/
      (1 + exp(T_AL/ Ti - T_AL/ T_L) + exp(T_AH/ T_H - T_AH/ Ti))
  }
  
  return(TC)
}
