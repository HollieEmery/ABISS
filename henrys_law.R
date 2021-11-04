
# Henry's law calculator ####
## written by Hollie Emery
## Updated 6/24/2021

# This function calculates the concentration (in mol/L) of dissolved gas in seawater or water,
# based on the (analytically determined or theoretical) concentration of gas at equilibrium with it,
# at given temperature and pressure conditions (default = ~lab conditions).

# For methane and hydrogen gas, a calculation that considers salinity is also available (and is the default). 

## Disclaimers and examples are given below (line ~120).

# Gases currently supported (*=supported with salinity adjustment): 
## hydrogen*, methane*, ethane, propane, butane, pentane, hexane,
## oxygen, nitrogen, nitrous oxide, argon, hydrogen sulfide, carbon dioxide. 


henry <- function(
  x,                 # headspace concentration in ppm (e.g. GC output)
  gas,               # analyte to calculate (see options above)
  TC = 22,           # temp in celcius (default = ~room temp)
  S = 34,            # salinity in PSU. Only used for salinity-adjusted gases (currently methane and hydrogen only)
  P = NA,            # pressure in atm. overrides any pressure calculation based on depth
  z = NA,            # depth in m, used to calculate P if P is not given
  temp_adj = TRUE,   # adjust for temperature 
  sal_adj = TRUE     # adjust for salinity if available [and temp] (currently methane and hydrogen only)
)
{
  # conversions
  TK = 273.15 + TC                 # celsius -> kelvin
  xp = x / 10^6                    # ppm -> proportion
  
  # set P if it's undefined
  if(is.na(P)){                    
    P = 1 +                        # 1 atm for the actual atmosphere
      ifelse(!is.na(z), z/10, 0)   # additional pressure from ocean depth @ 1 atm each 10 m
  }
  
  if(gas == "methane" & sal_adj == TRUE & temp_adj == TRUE){
    # coefficients and calculations from Wiesenburg & Guinasso 1979 (DOI 10.1021/je60083a006)
    A = c(-415.2807, 596.8104, 379.2599, -62.0757)  # nmol/L-atm
    B = c(-0.059160, 0.032174, -0.0048198) # nmol/L-atm
    
    # solubility constant
    KH = exp( A[1] +   A[2] * (100/TK) + A[3] * log(TK/100) + A[4] * (TK/100) +
                S*(B[1] + B[2] * (TK/100) + B[3] * (TK/100)^2) ) # nmol/(L-atm)
    KH = KH/10^9 #mol/(L-atm)
    
  }else if(gas == "hydrogen" & sal_adj == TRUE & temp_adj == TRUE){
    # coefficients and calculations from Wiesenburg & Guinasso 1979 (DOI 10.1021/je60083a006)
    A = c(-317.4669, 455.8526, 297.5313, -49.2778)  # nmol/L-atm
    B = c(-0.070143, 0.041069, -0.0063763) # nmol/L-atm
    
    # solubility constant
    KH = exp( A[1] +   A[2] * (100/TK) + A[3] * log(TK/100) + A[4] * (TK/100) +
                S*(B[1] + B[2] * (TK/100) + B[3] * (TK/100)^2) ) # nmol/(L-atm)
    KH = KH/10^9 #mol/(L-atm)
    
  }else{
    # Coefficient values are the best consensus (or most recent) literature values from Sander 2015 
    # and agree with NIST values.
    tab = data.frame( 
      Ko = c( 
        1.4e-5,                     # methane
        1.9e-5,                     # ethane
        1.5e-5,                     # propane
        1.2e-5,                     # butane
        8.0e-6,                     # pentane
        6.0e-6,                     # hexane
        7.8e-6,                     # hydrogen
        1.3e-5,                     # oxygen
        6.4e-6,                     # nitrogen
        2.4e-4,                     # nitrous oxide
        1.4e-5,                     # argon
        1.0e-3,                     # hydrogen sulfide
        3.3e-4                      #carbon dioxide NOTE: dissolved CO2 gas, NOT total DIC! You need pH for the whole picture.
      ) * 101.325,             ## mol/(m^3*Pa) -> mol/(L*atm) 
      
      dT = c(
        1900,                       # methane  
        2400,                       # ethane
        2700,                       # propane
        3100,                       # butane
        3400,                       # pentane
        3800,                       # hexane
        530,                        # hydrogen
        1500,                       # oxygen
        1300,                       # nitrogen
        2600,                       # nitrous oxide
        1500,                       # argon
        2100,                       # hydrogen sulfide
        2400                        # carbon dioxide 
        )
    ) 
    
    row.names(tab) = c("methane","ethane","propane","butane","pentane","hexane",
                       "hydrogen","oxygen","nitrogen","nitrous oxide","argon",
                       "hydrogen sulfide","carbon dioxide")
    
    # solubility constant with or without temperature correction
    if(temp_adj == TRUE){
      KH = tab[gas,"Ko"] * exp(tab[gas,"dT"] * (1/TK - 1/298.15))  # mol/(L-atm)
    }else{
      KH = tab[gas,"Ko"] # mol/(L-atm)
    }
    
  }
  
  # concentration
  Cw = (KH * xp * P )          # mol/L

  return(Cw)
  
}


# Examples ####
henry(2, "methane", S=0)                               # 2 ppm methane (~atmospheric levels) in DI water in the lab
henry(1000000, "methane", TC=4, z=1000)                # 100% methane on the seafloor - a cold seep perhaps?
henry(500000, "hydrogen", P=100)                       # 50% hydrogen @ 100 atm
henry(500000, "hydrogen", P=100, sal_adj=FALSE)        # suppress the salinity adjustment for comparison
henry(9340, "argon", TC=10)                            # argon (~atmospheric levels) near the sea surface
henry(.21*1e6, "oxygen", TC=16)                        # O2 (~atmospheric level of 21%) near sea surface
henry(.21*1e6, "oxygen", TC=16)*1000*1000              # convert to micromolar

# Disclaimers ####
# 1) In general this function is most accurate for values of temperature and pressure that are nearish to STP.
# 2) Salinity usually lowers gas solubility, so this function will given an overestimate for gases that don't have salinity adjustment.
# 3) Carbon dioxide is supported but use caution because this is dissolved CO2 gas, not the same as DIC. You need to know pH
#    to do the calculations that will tell you the whole picture of dissolved carbon dioxide because of bicarbonate buffering.

# References ####
# 1) Wiesenburg, DA & Guinasso, NL, 1979. Equilibrium solubilities of methane, carbon monoxide, and hydrogen in water and sea water
#    Journal of Chemical and Engineering Data 24(4) 356-360. DOI:10.1021/je60083a006
# 2) Sander, R, 2015. Compilation of Henry's law constants (version 4.0) for water as solvent.
#    Atmospheric Chemistry and Physics 16(8):4399-4981. DOI: 10.5194/acp-15-4399-2015. also available at www.henrys-law.org

