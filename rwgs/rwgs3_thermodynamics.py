from numpy import *
from scipy.optimize import *
import random
#from sympy.core.symbol import symbols
#from sympy.solvers.solveset import nonlinsolve

R = 8.3145 #J/mol/K
T = 600 + 273.15 #K
P = 1.01325*10**5 #J/m^3
n0 = {'CO2': 1, 'H2': 3, 'CO': 0, 'H2O': 0, 'CH4': 0} #kmol

#std Gibbs FE of formation
Gf = {'CO2': -394.38, 'H2': 0, 'CO': -137.16, 'H2O': -228.614, 'CH4': -50.45}
#rxn1: Reverse Water Gas Shift
G1 = Gf['CO'] + Gf['H2O'] - Gf['CO2'] - Gf['H2']
#rxn2: Methanation
G2 = Gf['CH4'] + Gf['H2O'] - Gf['CO'] - 3*Gf['H2']
#rxn3: the Sabatier rxn
G3 = Gf['CH4'] + 2*Gf['H2O'] - Gf['CO2'] - 4*Gf['H2']

#Calculating equilibrium constants K1, K2, and K3
K1 = exp(-G1/R/T)
K2 = exp(-G2/R/T)
K3 = exp(-G3/R/T)

def GasEq(zeta):
    #z = extent of rxn (kmol)
    z1 = zeta[0]
    z2 = zeta[1]
    z3 = zeta[2]
    print('This is the iterative zeta:' + str(zeta))
    #ntot = ((n0['H2']  - z1 - 3*z2 - 4*z3) + (n0['CO2'] - z1 - z3) + (n0['CO']  + z1 - z2) + (n['H2O'] = n0['H2O'] + z1 + z2 + 2*z3) + (n0['CH4'] + z2 + z3))

    #system of equations, equilibrium constants in terms of partial pressures
    E = [0]*3
    #E[0] = rwgs.
    E[0] = n0['CO']  + z1 - z2
    E[0] = E[0]*(n0['H2O'] + z1 + z2 + 2*z3)
    E[0] = E[0]/(n0['CO2'] - z1 - z3)
    E[0] = E[0]/(n0['H2']  - z1 - 3*z2 - 4*z3)
    E[0] = E[0] - K1
    #E[1] = methanation
    E[1] = (1/(P**2)*((n0['H2']  - z1 - 3*z2 - 4*z3) + (n0['CO2'] - z1 - z3) + (n0['CO']  + z1 - z2) + (n0['H2O'] + z1 + z2 + 2*z3) + (n0['CH4'] + z2 + z3))**2)
    E[1] = E[1]*(n0['CH4'] + z2 + z3)
    E[1] = E[1]*(n0['H2O'] + z1 + z2 + 2*z3)
    E[1] = E[1]/(n0['CO']  + z1 - z2)
    E[1] = E[1]/((n0['H2']  - z1 - 3*z2 - 4*z3)**3)
    E[1] = E[1] - K2
    #E[2] = Sabatier reaction
    E[2] = (1/(P**2)*((n0['H2']  - z1 - 3*z2 - 4*z3) + (n0['CO2'] - z1 - z3) + (n0['CO']  + z1 - z2) + (n0['H2O'] + z1 + z2 + 2*z3) + (n0['CH4'] + z2 + z3))**2)
    E[2] = E[2]*(n0['CH4'] + z2 + z3)
    E[2] = E[2]*((n0['H2O'] + z1 + z2 + 2*z3)**2)
    E[2] = E[2]/(n0['CO2'] - z1 - z3)
    E[2] = E[2]/((n0['H2']  - z1 - 3*z2 - 4*z3)**4)
    E[2] = E[2] - K3
    print('This is the iterative E:' + str(E))
    return E

zetaGuess = array([random.uniform(.1,.5), random.uniform(0,.1), random.uniform(0,.2)])

zeta = fsolve(GasEq, zetaGuess)
print('This is the final zeta:' + str(zeta))

n = {}
n['H2']  = n0['H2']  - zeta[0] - 3*zeta[1] - 4*zeta[2]
n['CO2'] = n0['CO2'] - zeta[0] - zeta[2]
n['CO']  = n0['CO']  + zeta[0] - zeta[1]
n['H2O'] = n0['H2O'] + zeta[0] + zeta[1] + 2*zeta[2]
n['CH4'] = n0['CH4'] + zeta[1] + zeta[2]
print('This is the final n:' + str(n))

#calculating total number of moles, and mole fractions of each component.
ntot = 0
for key in n:
    ntot = ntot + n[key]
y = {}
for key in n:
    y[key] = n[key] / ntot
print('This is the final y:' + str(y))


#plt.plot(x=T, y['CO2'], color = "yellow")
#plt.plot(x=T, y['H2'], color = "green")
#plt.plot(x=T, y['CO'], color = "violet")
#plt.plot(x=T, y['H2O'], color = "blue")
#plt.plot(x=T, y['CH4'], color = "red")
#plt.show()
