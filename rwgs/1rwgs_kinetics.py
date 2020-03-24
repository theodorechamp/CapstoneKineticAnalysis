import matplotlib.pyplot as plt
import numpy as numpy
import math
from scipy.integrate import solve_ivp

#molar flow rates of each species (for Ergun equation)
#def nflow(p, dV, R, T):
#    nflow = [0]*len(p)
#    for i in range(len(p)):
#        nflow[i] = p[i]*dV/R/T
#    return(nflow)

#total mass flow rate (for mass flux in Ergun equation)
def mflow(n, MW):
    mflowTot = 0; #
    for i in range(len(MW)):
        mflowTot = mflowTot + n[i]*MW[i]
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

# f(t, y, c)
# f(V, y, P0CO2, P0H2, k1L0, K2, K3, K, X, p, dV, MW, T, R, mflowTot, D)
# y = [nCO2, nH2, nCO, nH2O, P]
def reaction(V, y, C):
    # unpacking C = [P0CO2, P0H2, k1L0, K2, K3, dV, MW, T, R, mflowTot, D]
    P0CO2, P0H2, k1L0, K2, K3, dV, MW, T, R, mflowTot, D = C
    #calculate p from y = [nCO2, nH2, nCO, nH2O, P]
    p = [0]*4
    for i in range(0,3):
        p[i] = y[i]*R*T/dV
    K = p[2]*p[3]/p[0]/p[1]
    X = (P0CO2 - p[0])/P0CO2
    #rate units = mol/gcat/hr
    r = k1L0*P0CO2*(P0H2*((1-X)**2) - P0CO2*(X**2)/K)
    r = r/(P0H2*(1-X) + (K2**0.5)*P0H2*1.5*(((1-X)**1.5)) + P0CO2*X/K2/K3)
    gCat = 85
    ntot = 0
    for i in range(0,3):
        ntot = ntot + y[i]
    Ptot = y[4]
    #first four entries are dn/dV, 5th is dPtot/dV
    dydV = [-1*r*gCat*Ptot/(ntot*R*T), -1*r*gCat*Ptot/(ntot*R*T),
            1*r*gCat*Ptot/(ntot*R*T), 1*r*gCat*Ptot/(ntot*R*T),
            ergun(p, MW, T, R, mflowTot, D)]
    return(dydV) #unsure of whether to return X, since we are going to solve for it


# OVERALL REACTION: CO2 + H2 -> CO + H2O
# We have rate r, [r]=mol/hr/g catalyst
# dnCO2/dt = -r*gCat
# dnH2/dt = -r*gCat
# dnCO/dt = r*gCat
# dnH2O/dt = r*gCat
#
# dn/dV = dn/dt/Vdot = dn/dt/ndot/R/T*Ptot = nu*r*gCat*Ptot/(ndot*R*T)
#
# y(t) = [nCO2[x1,x2,x3], nH2, nCO, nH2O]
#


def main():
    #constant parameters
    R = 0.08205 #gas constant. units = L atm K^-1 mol^-1
    T = 523 #temperature. units = K
    #molecular weight dictionary for all species
    MWDict = {'CO2': 44.008, 'H2': 1.008, 'CO': 28.008, 'H2O': 18.016}
    MW = list(MWDict.values()) #create list instead of dictionary for calculations

    #dictionary for feed molar flowrates of each species
    nDict = {'CO2':0.25, 'H2': 0.74998, 'CO': 0.00001, 'H2O': 0.00001}
    n = list(MWDict.values())
    ntot = 0
    for i in range(len(n)):
        ntot = ntot + n[i]

    mflowTot = mflow(n, MW)
    #feed mole fraction dictionary.
    yDict = {'CO2': nDict['CO2']/ntot, 'H2': nDict['H2']/ntot, 'CO': nDict['CO']/ntot, 'H2O': nDict['CO2']/ntot}
    Ptot = 33*1.01325 #feed pressure. units converted from bar to atm
    #initial partial pressures. units = atm
    pDict = {'CO2': Ptot*yDict['CO2'],
            'H2': Ptot*yDict['H2'],
            'CO': Ptot*yDict['CO'],
            'H2O': Ptot*yDict['H2O']
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
#    dV = V[1]
#    print('this is dV: ' + str(dV))
#    flow = nflow(p, dV, R, T)
#    print('CO2 nflow is ' + str(flow[0]))
#    mflowTot = mflow(flow, MW)
#    print('mflowTot is ' + str(mflowTot))
#    dPdV = ergun(p, MW, T, R, mflowTot, D)
#    print('dPdV is ' + str(dPdV))
#    print(p)

    # y0 = [nCO20, nH20, nCO0, nH2O0, P0]
    C = [P0CO2, P0H2, k1L0, K2, K3, V[1], MW, T, R, mflowTot, D]
    y0 = [ nDict['CO2'], nDict['H2O'], nDict['CO'], nDict['H2O'], Ptot]
    #(V, y, P0CO2, P0H2, k1L0, K2, K3, dV, MW, T, R, mflowTot, D)
    sol = solve_ivp(lambda t, y: reaction(t,y,C), [V[0], V[-1]], y0, t_eval = V)


if __name__ == '__main__':
    main()
