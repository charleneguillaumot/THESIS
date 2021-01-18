#'@export
beta0 <- function (x0,x1){
  x03 = x0^(1/3); x13 = x1^(1/3); a3 = sqrt(3)
  f1 = -3 * x13 + a3 * atan((1+2 * x13)/a3) - log (as.complex(x13 - 1))+
    log(1 + x13 + x13^2)/2
  f0 = -3 * x03 + a3 * atan((1+2 * x03)/a3) - log (as.complex(x03 - 1))+
    log(1 + x03 + x03^2)/2
  f = f1 - f0
  return(f)
} 

#'@export
dget_l_ISO <- function (vH, l, parms) {  # growth for isometric adult)
  if (length(parms)< 3){
    warning ("parms vector must include at least k,g,f (same order)")
    stop()
  }
  if (length(parms)==3){
    k = parms[1]
    g = parms[2]
    f = parms[3]
    sM = 1
    }else{
      k = parms[1]
      g = parms[2]
      f = parms[3]
      sM = parms[4]
      }

  r = g * (f * sM - l)/ l / (f + g)
  dl = l * r / 3
  dvH = f * l^2 * (sM - l * r / g) - k * vH
  dl = dl / dvH
  return(list(dl))
}  # get length asuming isomorphic growth (VB)

#'@export
dget_l_V1 <- function (vH, l, parms){   # growth for v1 morph (for length at metamorphosis)
  k = parms[1]
  g= parms[2]
  f=parms[3]
  lref=parms[4]
  r=((g*(f - lref))/lref)/(g+f) # specific growth rate
  dl = r * l/3 # d/dt l 
  dvH = f * l^3 * (1/lref - r/g) - k * vH; # d/dt vH
  dl = dl/dvH;
  return(list(dl))
}   # get length at metamorphosis asumes V1-morph growth

#'@export
dget_LUH <- function (a, LUH, p) {
  kap = p$kap
  v = p$v
  kJ= p$k_J
  g = p$g
  Lm = p$L_m
  L = LUH[1]  # cm, structural length
  U = LUH[2]  # d cm^2, scaled reserves M_E/{J_EAm}
  H = LUH[3]  # d cm^2, scaled maturity M_H/{J_EAm}
  eL3 = U * v
  gL3 = g * L^3
  SC = L^2 * (1 + L/(g*Lm))*g*eL3/(gL3 + eL3)
  dL = v * (eL3 - L^4 / Lm)/(3 * (eL3+gL3))
  dU = -SC
  dH = (1-kap)* SC - kJ * H
  return(list(c(dL,dU,dH)))
}  # get structural length and reserves at age during embryonic development

#'@export
get_uE0 <- function (p, e){ # initial scaled reserves
  e=e        # -, scaled functional response for the mother
  g = p$g    # -, energy investment ratio
  v = p$v    # cm2/d, energy conductance
  k_M=p$k_M  # 1/d, somatic maintenance rate coefficient
  k = p$k    # -, somatic/reproductive ratio
  v_Hb = p$v_Hb 
  
  ##### calculating lb --> scaled length at birth ####
  
  n = 1000 + round(1000 * max(0, k - 1))
  xb = g/ (g + e)
  xb3 = xb ^ (1/3)
  x = seq(1e-6, xb, length=n)
  dx = xb/ n
  x3 = x ^ (1/3)
  lb = v_Hb^(1/3) # exact solution for k==1, to be used as starting value.
  b = beta0(x, xb)/ (3 * g)
  t0 = xb * g * v_Hb
  i = 0
  norm = 1 # make sure that we start Newton Raphson procedure
  ni = 100 # max number of iterations
  
  while ((i < ni) & (norm > 1e-8)) {
    l = x3 / (xb3/ lb - b);
    s = (k - x) / (1 - x) * l/ g / x;
    y = exp( - dx * cumsum(s))
    vb = y[n];
    r = (g + l)
    rv = r / y;
    t = t0/ lb^3/ vb - dx * sum(rv);
    dl = xb3/ lb^2 * l^2 / x3
    dlnl = dl / l;
    dv = y * exp( - dx * cumsum(s * dlnl));
    dvb = dv[n]
    dlnv = dv / y
    dlnvb = dlnv[n];
    dr = dl 
    dlnr = dr / r;
    dt = - t0/ lb^3/ vb * (3/ lb + dlnvb) - dx * sum((dlnr - dlnv) * rv);
    lb = lb - t/ dt # Newton Raphson step
    lb = Re(lb)
    norm = Re(t^2)
    i = i + 1   
  }
  
  xb = g/ (e + g)
  u0 = Re( (3 * g/ (3 * g * xb^(1/ 3)/ lb - beta0(0, xb)))^3)
  uE0 = u0 * v^2/ g^2/ k_M^3 # d.cm^2, initial scaled reserve
  return(uE0)
}

#'@export
random.sample<-function(n,mean,sd,lowerBound,upperBound){
  range <- upperBound - lowerBound
  m <- (mean-lowerBound) / range #mapping mean to 0-1 range
  s <- sd / range #mapping sd to 0-1 range
  a <- (m^2 - m^3 - m*s^2)/s^2 #calculating alpha for rbeta 
  b <- (m-2*m^2+m^3-s^2+m*s^2)/s^2 #calculating beta for rbeta
  data <- rbeta(n,a,b)  #generating data
  data <- lowerBound + data * range #remaping to given bounds
  return(data)
}
