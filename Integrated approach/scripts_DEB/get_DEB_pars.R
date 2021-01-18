#'@export
#'@importFrom "R.matlab" "readMat"

get_DEB_pars<- function (file) {
  # get pars. This function reads the MatLab.mat file and gives a list with all
  # the DEB parameters.Parameters can also be read from a csv file filled as the
  # available template (R-package data).the only input needed is the name of the
  # species as "Genus_species" format (make sure it is a string). Those
  # parameters are fixed and do not depend on food levels and they are given at
  # a reference temperature (T_ref) included in the T.pars. the function has
  # different outputs: p = a list containing all the DEB parameters also
  # contains a vector t.pars containing the Arrhenius temperature parameters for
  # temperature correction (T_ref,T_A,T_L,T_H,T_AL,T_AH), according to
  # disponibility and also the matrices eta.O and n.M for the coupling of energy
  # and mass.
#   T_ref     K           Temperature for which parameters are given
#   z         -           zoom factor; for z=1: L_m = 1cm
#   F_m       l/d.cm^2
#   kap_X     -           digestion efficiency of food to reserves
#   kap_P     -           defacation efficiency
#   v         cm/d        energy conductance
#   kap       -           allocation fraction to soma = growth + somatic maintenance
#   kap_R     -           reproduction efficiency
#   p_M       J/d.cm^3    vol-specific somatic maintenance
#   p_T       J/d.cm^2    surface-specific somatic maintenance
#   k_J       1/d         maturity maintenance rate coefficien
#   E_G       J/cm^3      specific cost for structure
#   E_Hb      J           Maturity at birth (start of exogenous feeding)
#   E_Hj      J           Maturity at metamorphosis
#   E_Hp      J           Maturity at puberty (start of allocation to Er: reproduction buffer)
#   h_a       1/d^2       Weibull aging acceleration
#   s_G       -           Gompertz stress coefficient
#   T_A       K           Arrhenius Temperature
#   T_L       K           Lower temperature limit
#   T_H       K           Upper temperature limit
#   T_AL      K           Arrhenius Temperature at lower limit
#   T_AH      K           Arrhenius Temperature at upper limit
#   del_M     -           post metamorphic shape coefficient
#   del_Mj    -           pre metamorphic shape coefficient

# check file extension - MatLab files are from the matlab DEB toolbox, txt/csv files are from this package template.
ext <- substring(file, nchar(file)-2)

if (ext == 'mat'){
  x=readMat(file)
  l=x$par
  l=l[,,1]
  n.par=names(l)
  f.par=gsub(".", "_", n.par, fixed=TRUE)
  p=list()
  # extract primary parameters from .mat file
  suppressWarnings(
    for (i in 1:length(n.par)) {
      p[f.par[i]]=unlist(l[n.par[i]])
    }
  )
}else{
  if (ext == 'csv'){
    x = read.table(file, header=TRUE, sep=';')
    x = as.list(x[,1:2])
    l = x$Value
    n.par = as.character(x$names)
    p=list()
    # extract primary parameters from .mat file
      for (i in 1:length(n.par)) {
        p[n.par[i]]=l[i]
      }
    # adding chemical and elemental parameters (they can be modified on system file)
    data("elemental_composition")
    chems = as.list(elemental_composition[,1:2])
    l = chems$Value
    n.par = as.character(chems$names)
    # extract primary parameters from .mat file
    for (i in 1:length(n.par)) {
      p[n.par[i]]=l[i]
    }
  }else{stop('wrong file type')}

}

# calculating aditional of DEB parameters
attach(p)

p["w_X"] = w_X = 23.9;
p["w_V"] = w_V = 23.9;
p["w_E"] = w_E = 23.9;               # g/mol,  molecular weights for food, structure and reserve

p["p_Am"] = p_Am = z * p_M/ kap;                # J/d cm^2, {p_Am} max spec assimilation flux
p["J_E_Am"] = J_E_Am = p_Am/ mu_E;              # mol/d cm^2, {J_EAm}, max surface-spec assimilation flux
p["k_M"] = k_M = p_M/ E_G;                      # 1/d, somatic maintenance rate coefficient
p["k"]= k = k_J/ k_M;                           # -, maintenance ratio

p["M_V"] = M_V = d_V/ w_V;                     # mol/cm^3, [M_V], volume-specific mass of structure
p["kap_G"] = kap_G = M_V * mu_V/ E_G;          # -, growth efficiency
p["y_V_E"] = y_V_E = mu_E * M_V/ E_G;          # mol/mol, yield of structure on reserve
p["y_E_V"] = y_E_V = 1/ y_V_E;                 # mol/mol, yield of reserve on structure
p["y_E_X"] = y_E_X = kap_X * mu_X/ mu_E;       # mol/mol, yield of reserve on food
p["y_X_E"] = y_X_E = 1/ y_E_X;                 # mol/mol, yield of food on reserve
p["y_P_X"] = y_P_X = kap_P * mu_X/ mu_P;       # mol/mol, yield of faeces on food
p["y_P_E"] = y_P_E = y_P_X/ y_E_X;             # mol/mol, yield of faeces on reserve

p["E_m"] = E_m = p_Am/ v;                      # J/cm^3, [E_m] reserve capacity
p["m_Em"] = m_Em = y_E_V * E_m/ E_G;           # mol/mol, reserve capacity
p["g"] = g = E_G / kap / E_m                   # -, energy investment ratio
p["w"] = w = m_Em * w_E/ w_V;                  # -, contribution of reserve to weight

p["L_m"] = L_m = v/ k_M/ g;                    # cm, maximum structural length
p["L_T"] = L_T = p_T/ p_M;                     # cm, heating length (also applies to osmotic work)
p["l_T"] = l_T = L_T/ L_m;                     # -, scaled heating length


# calculating values of maturities

if (exists('E_Hb')){
  # maturity at birth
  p["M_Hb"] = M_Hb = E_Hb/ mu_E;               # mol, maturity at birth
  p["U_Hb"] = U_Hb = M_Hb/ J_E_Am;             # cm^2 d, scaled maturity at birth
  p["u_Hb"] = u_Hb = U_Hb * g^2 * k_M^3/ v^2;  # -, scaled maturity at birth
  p["V_Hb"] = V_Hb = U_Hb/ (1 - kap);          # cm^2 d, scaled maturity at birth
  p["v_Hb"] = v_Hb = V_Hb * g^2 * k_M^3/ v^2;  # -, scaled maturity at birth
}

if (exists('E_Hj')){  # abj model!
  # maturity at metamorphosis
  p["M_Hj"] = M_Hj = E_Hj/ mu_E;               # mol, maturity at metamorphosis
  p["U_Hj"] = U_Hj = M_Hj/ J_E_Am;             # cm^2 d, scaled maturity at metamorphosis
  p["u_Hj"] = u_Hj = U_Hj * g^2 * k_M^3/ v^2;  # -, scaled maturity at metamorphosis
  p["V_Hj"] = V_Hj = U_Hj/ (1 - kap);          # cm^2 d, scaled maturity at metamorphosis
  p["v_Hj"] = v_Hj = V_Hj * g^2 * k_M^3/ v^2;  # -, scaled maturity at metamorphosis
}

if (exists('E_Hp')){
  # maturity at puberty
  p["M_Hp"] = M_Hp = E_Hp/ mu_E;               # mol, maturity at puberty
  p["U_Hp"] = U_Hp = M_Hp/ J_E_Am;             # cm^2 d, scaled maturity at puberty
  p["u_Hp"] = u_Hp = U_Hp * g^2 * k_M^3/ v^2;  # -, scaled maturity at puberty
  p["V_Hp"] = V_Hp = U_Hp/ (1 - kap);          # cm^2 d, scaled maturity at puberty
  p["v_Hp"] = v_Hp = V_Hp * g^2 * k_M^3/ v^2;  # -, scaled maturity at puberty
}

#mass-power couplers
p["eta_XA"] = eta_XA = y_X_E/mu_E;             # mol/J, food-assim energy coupler
p["eta_PA"] = eta_PA = y_P_E/mu_E;             # mol/J, faeces-assim energy coupler
p["eta_VG"] = eta_VG = y_V_E/mu_E;             # mol/J, struct-growth energy coupler

n_O=matrix(c(p$n_CX,p$n_CV,p$n_CE,p$n_CP,
             p$n_HX,p$n_HV,p$n_HE,p$n_HP,
             p$n_OX,p$n_OV,p$n_OE,p$n_OP,
             p$n_NX,p$n_NV,p$n_NE,p$n_NP)
           , nrow=4, ncol=4, byrow=TRUE)
p["n_O"]=list(n_O)

eta_O = matrix(                 # mol/J, mass-energy coupler
  c(-eta_XA,0,0,0,0, eta_VG, 1/mu_E, -1/mu_E, -1/mu_E, eta_PA,0,0),
  nrow=4,
  ncol=3,
  byrow=TRUE)
p["eta_O"]=list(eta_O)

n_M = matrix(
  c(1,0,0,0,0,2,0,3,2,1,2,0,0,0,0,1),
  nrow=4,
  ncol=4,
  byrow=TRUE)
p["n_M"]=list(n_M)

# Arrhenius temperature parameters for temperature corrections
t.pars=c()
t.pars[1]=T_ref
t.pars[2]=T_A

if (exists("T_AH")){
  t.pars[3]=T_L
  t.pars[4]=T_H
  t.pars[5]=T_AL
  t.pars[6]=T_AH
} else {
  if (exists("T_AL")){
    t.pars[3]=T_L
    t.pars[4]=T_AL
  }
}

p["t.pars"]=list(t.pars)

detach(p)

return (p)
}
