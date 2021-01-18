#'@export
get_powers_j <- function(p, L, e, s_M=0){
  ## p - a vector containing the DEB parameters
  ## L - structural length of the individual to calculate dER
  ## e - scaled reserves (=f when in equilibrium with available food)
  ## s_M, optional acceleration factor.
  
  ## This function calculates scaled powers assimilation, mobilisation, somatic maintenance, maturity maintenance,
  ## growth, reproduction and dissipation as function of length and scaled reserves
  
  #Output
  #
  # pACSJGRD: (n,7)- vector or matrix with scaled powers p_A, p_C, p_S, p_J, p_G, p_R, p_D / L_m^2 {p_Am}
  
  
  if (s_M == 0) {
    tj= get_tj(p, e)
    s_M=tj$lj/tj$lb
  }
 
  
  l = L/p$L_m
  pA = e * s_M * l^2;                         # assim
  pC = l^2 * ((p$g + p$l_T) * s_M + l )/ (1 + p$g/ e);   # mobilisation
  pS = p$kap * l^2 * (l + p$l_T);               # somatic  maint
  pJ = p$k * p$u_Hp;                           # maturity maint
  pG = p$kap * pC - pS;                        # growth
  pR = (1 - p$kap) * pC - pJ;                  # maturation/reproduction
  pD = pS + pJ + (1 - p$kap_R) * pR ;          # dissipation
  
  power = cbind(pA, pC, pS, pJ, pG , pR, pD);
  return(power)
}
