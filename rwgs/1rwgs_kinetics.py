import matplotlib.pyplot as plt
import numpy as numpy
import math
from scipy.integrate import solve_ivp

#molar flow rates of each species (for Ergun equation)
def nflow(p, dV, R, T):
    nflow = [0]*len(p)
    for i in range(len(p)):
        nflow[i] = p[i]*dV/R/T
    return(nflow)
#total mass flow rate (for mass flux in Ergun equation)
def mflow(nflow, MW):
    mflowTot = 0; #
    for i in range(len(MW)):
        mflowTot = mflowTot + nflow[i]*MW[i]
    return(mflowTot)
#Ergun equation
def ergun(p, MW, T, R, mflowTot, D):
    rho = 0 #total density of the gas mixture
    for i in range(len(p)):
        rho = rho + p[i]*MW[i]/(T*R)
    phi = 0.3 #void fraction
    Dp = D/8 #Particle diameter. Using heuristic Weimer gave in class on 10/10/19 that particle diameter should be 1/8 of tube diameter
    Ac = math.pi*(D**2)/4 #cross sectional area
    G0 = mflowTot/Ac #mass flux
    mu = 3.6*0.01849 #3.6 = unit conversion. Value from HYSYS. Units = kg m^-1 h^-1
    #differential pressure change
    dPdV = -150/Ac*(mu*(1-phi)/(Dp*G0)+(7/4))*((1-phi)/phi**3)*(G0**2/(rho*Dp))
    return(dPdV)

#def reaction(p, k1L0, K2, K3, K, X)
    # rate = k1L0*p0CO2*(p0H2*((1-X)**2) - p0CO2*(X**2)/K)
    # rate = rate/(p0H2*(1-X) + (K2**0.5)*p0H2*1.5*(((1-X)**1.5)) + p0CO2*X/K2/K3)
    # p[0] = p[0] - rate
    # K = p['CO']*p['H2O']/p['CO2']/p['H2']
    #return(rate, X) #unsure of whether to return X, since we are going to solve for it

def main():
    #constant parameters
    R = 0.08205 #gas constant. units = L atm K^-1 mol^-1
    T = 523 #temperatire. units = K
    #molecular weight dictionary for all species
    MWDict = {'CO2': 44.008, 'H2': 1.008, 'CO': 28.008, 'H2O': 18.016, 'CH4': 16.04}
    MW = list(MWDict.values()) #create list instead of dictionary for calculations

    #feed mole fraction dictionary.
    yDict = {'CO2': 0.25, 'H2': 0.74997, 'CO': 0.00001, 'H2O': 0.00001, 'CH4': 0.00001}
    P = 33*1.01325 #feed pressure. units converted from bar to atm
    #initial partial pressures. units = atm
    pDict = {'CO2': P*yDict['CO2'],
            'H2': P*yDict['H2'],
            'CO': P*yDict['CO'],
            'H2O': P*yDict['H2O'],
            'CH4': P*yDict['CH4']
            }
    #extract initial CO2 and H2 pressures (needed for rate calculations).
    P0CO2 = pDict['CO2']
    P0H2 = pDict['H2']
    p = list(pDict.values()) #create list instead of dictionary for calculations

    #kinetic parameters. taken from results of Gines et al., 1997
    #k1L0 = (rate of fwd CO2 dissoc.)*(conc. of empty catalyst sites)
    k1L0 = 0.1287 #units = mol/h/gcat/atm**2
    K2 = 0.039 #equilibrium constant for H2 adsoprtion
    K3 = 0.00000125 #equilibrium constant for H2O dissociation
    X = 0.0000001 #conversion (initial)
    K = pDict['CO']*pDict['H2O']/pDict['CO2']/pDict['H2'] #overall equilibrium constant
    rate = k1L0*pDict['CO2']/(1 + (K2*pDict['H2'])**0.5) #RWGS rxn rate

    #reactor dimensions.
    D = 0.0254 #tube diameter. units = m
    L = 1 #tube length. units = m
    Vtot = 1000*L*math.pi*(D**2)/4 #total volume of 1 tube. units converted to L
    V = numpy.linspace(0,Vtot,101)
    #for now, set dV as a small value. will later have dV be iteration through V array
    dV = V[1]
    print('this is dV: ' + str(dV))
    flow = nflow(p, dV, R, T)
    print('CO2 nflow is ' + str(flow[0]))
    mflowTot = mflow(flow, MW)
    print('mflowTot is ' + str(mflowTot))
    dPdV = ergun(p, MW, T, R, mflowTot, D)
    print('dPdV is ' + str(dPdV))
    print(p)


if __name__ == '__main__':
    main()
