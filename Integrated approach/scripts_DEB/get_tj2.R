#'@import pracma
#'@import deSolve

#'@export
get_tj <- function (p, f) {
  
  #### this script calculate some of the variables associated to the life cycle of
  #### the organism is a equivalent to get_tj but restricted to the particular
  #### case of abj.
  
  
  ##### Obtaining lb, lj, lp, li, rj, rB, s_M, and tb, tj, tp, ti (ti is the age
  ##### it takes to almost reach li, not to be mistaken by the tm (mean life) of
  ##### the DEB model).
  
  # input --> p: a list containing all the DEB parameters
  # f --> scaled functional response for which these parameters should be calculated
  
  require(pracma)
  require(deSolve)
  
  if (f>1 || f<=0){
    warning("f must be a value between 0 and 1")
    stop()
  }
  
  if (length(p)< 4){
    warning ("p vector must include at least g, k, v_Hb and v_Hj (same order)")
    stop()
  }
  
  if (length(p) > 5) { attach(p)     # attaching the names of the parameters contained in p for their use in this function
  } else {
    g = p[1]
    k = p[2]
    v_Hb = p[3]
    v_Hj = p[4]
    v_Hp = 0
    if (length(p)==5){
      v_Hp = p[5]
    }
  }  
    

  ##### calculating lb --> scaled length at birth ####
  
  
  n = 1000 + round(1000 * max(0, k - 1))
  xb = g/ (g + f)
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
  
  ##### obtaining tb --> scaled age at birth ####
  
  dget_tb <- function(x, ab, xb){
    f = x^(-2/3) / (1 - x) / (ab - beta0(x, xb))
    return(f)
  }  # get scaled age at birth
  
  xb = g/ (f + g);
  ab = 3 * g * xb^(1/ 3)/ lb
  tb = 3 * quad(dget_tb, 1e-15, xb,  tol=1e-15, trace = FALSE, ab, xb);
  tb = Re(tb)
  
  
  ##### calculating lj and s_M --> scaled length at metamorphosis and acceleration factor ####
    
  vH = seq(v_Hb,v_Hj, length.out=45)
  vH_lj = ode(lb,vH,dget_l_V1, c(k,g,f,lb), method='ode45')
  lj= vH_lj[45,2] 
  s_M = lj/lb # -, acceleration factor
  
  
  ##### calculating scaled ultimate length li; scaled ages tj and growth rates rj and rB  #####
  # li scaled ultimate length
  # rj scaled exponential growth rate between b and j
  # tj scaled age at metamorphosis
  # rB scaled von Bert growth rate between j and i
  # ti scaled age at a value close to li.
  
  li = f * s_M  
  rj = g * (f/ lb - 1 - l_T/ lb)/ (f + g) 
  tj = tb + log(s_M) * 3/ rj;            
  rB = 1/ 3/ (1 + f/ g);  
  ti = tj + log((li - lj)/ (li - li*0.99))/ rB
  
  ## create output
  
  ltr = list()  #named list containign the sought values (lb,lj,...)
  ltr["lb"]=lb   # scaled length at birth
  ltr["lj"]=lj   # scaled length at metamorphosis
  ltr["li"]=li   # scaled ultimate length (vB)
  ltr["tb"]=tb   # scaled age at birth
  ltr["tj"]=tj   # scaled age at metamorphosis
  ltr["ti"]=ti   # scaled age at which animal reach li
  ltr["s_M"]=s_M # acceleration factor
  ltr["rj"]=rj   # scaled growth rate at exponential phase
  ltr["rB"]=rB   # scaled von Bertalanffy growth rate
  
  
  # ##### calculate lp and tp --> scaled length and age at puberty #####
  # if (v_Hp > 0){ # corresponds to the matlab code at "if empty, no values are indicated for lp and tp
  # vH = seq (v_Hj, v_Hp, length.out=45)
  # vH_lp = ode( lj, vH, dget_l_ISO, c(k,g,f,s_M), method='ode45')
  # lp = vH_lp[45,2]
  # ltr["lp"]=lp   # scaled length at puberty
  # tp = tj + log((li - lj)/ (li - lp))/ rB
  # ltr["tp"]=tp   # scaled age at puberty
  # 
  # }
  
  detach(p)

  return(ltr)
}

