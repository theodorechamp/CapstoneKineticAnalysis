import numpy as numpy
from scipy.optimize import *
import math

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
def ergun(p, MW, T, R, mflowTot):
    rho = 0 #total density of the gas mixture
    for i in range(len(p)):
        rho = rho + p[i]*MW[i]/(T*R)
    phi = 0.3 #void fraction
    D = 0.0254 #tube diameter
    Dp = D/8 #Particle diameter. Using heuristic Weimer gave in class on 10/10/19 that particle diameter should be 1/8 of tube diameter
    Ac = math.pi*(D**2)/4 #cross sectional area
    G0 = mflowTot/Ac #mass flux
    mu = 3.6*0.01849 #3.6 = unit conversion. Value from HYSYS. Units = kg m^-1 h^-1
    #differential pressure change
    dPdV = -150/Ac*(mu*(1-phi)/(Dp*G0)+(7/4))*((1-phi)/phi**3)*(G0**2/(rho*Dp))
    return(dPdV)


#def main():
#constant parameters
R = 8.3144621 #L kPa K^-1 mol^-1
T = 673 #K
#molecular weight dictionary for all species
MWDict = {'CO2': 44.008, 'H2': 1.008, 'CO': 28.008, 'H2O': 18.016, 'CH4': 16.04}
MW = list(MWDict.values())

#partial pressures (p0 = initial, p = current)
p0Dict = {'CO2': 1, 'H2': 1, 'CO': 0.001, 'H2O': 0.001, 'CH4': 0.001}
pDict = p0Dict
p = list(pDict.values())

#kinetic parameters
k1 = 1          #rate of fwd CO2 dissociation
L0 = 1          #cocentration of empty catalyst sites
X  = 0.0000001     #conversion
K2 = 1          #equilibrium constant for H2 adsoprtion
K3 = 1          #equilibrium constant for H2O dissociation
K = pDict['CO']*pDict['H2O']/pDict['CO2']/pDict['H2']                               #overall equilibrium constant
rate = k1*L0*p0Dict['CO2']/(1 + (K2**0.5)*p0Dict['H2'])                     #RWGS rxn rate


D = 0.0254
L = 1
Vtot = L*math.pi*(D**2)/4
V = numpy.linspace(0,Vtot,101)
#for now, set dV as a small value. will later have dV be iteration through V array
dV = V[1]
print('this is dV: ' + str(dV))
nflow = nflow(p, dV, R, T)
print('CO2 nflow is ' + str(nflow[0]))
mflowTot = mflow(nflow, MW)
print('mflowTot is ' + str(mflowTot))
dPdV = ergun(p, MW, T, R, mflowTot)
print('dPdV is ' + str(dPdV))

    #include in iteration:
    #rxn rate, in terms of initial partial pressures of CO2 and H2. variable = X
    # rate = k1*L0*p0['CO2']*(p0['H2']*((1-X)**2) - p0['CO2']*(X**2)/K)
    # rate = rate/(p0['H2']*(1-X) + (K2**0.5)*p0['H2']*1.5*(((1-X)**1.5)) + p0['CO2']*X/K2/K3)
    #
    # K = p['CO']*p['H2O']/p['CO2']/p['H2']
    #
    # print('This is the rate: ' + str(rate))
    # print('This is K: ' + str(K))


# if __name__ == "__main__":
#     main()
