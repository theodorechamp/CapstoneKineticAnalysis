import math
from sympy.core.symbol import symbols
from sympy.solvers.solveset import nonlinsolve


R = 8.3145 #J/mol/K
T = 298.15 #K
P = 1.01325*10**5 #J/m^3

#Calculating equilibrium constants K1, K2, and K3
Gf = {'CO2': -394.38, 'H2': 0, 'CO': -137.16, 'H2O': -228.614, 'CH4': -50.45}
#rxn1: Reverse Water Gas Shift
G1 = Gf['CO'] + Gf['H2O'] - Gf['CO2'] - Gf['H2']
#rxn2: Methanation
G2 = Gf['CH4'] + Gf['H2O'] - Gf['CO'] - 3*Gf['H2']
#rxn3: the Sabatier rxn
G3 = Gf['CH4'] + 2*Gf['H2O'] - Gf['CO2'] - 4*Gf['H2']

K1 = math.exp(-G1/R/T)
K2 = math.exp(-G2/R/T)
K3 = math.exp(-G3/R/T)

print(K1)
print(K2)
print(K3)
