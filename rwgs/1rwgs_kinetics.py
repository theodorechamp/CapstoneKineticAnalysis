from numpy import *
from scipy.optimize import *


#initial parameters
k1 = 1          #rate of fwd CO2 dissociation
L0 = 1          #cocentration of empty catalyst sites
X  = 0.0001     #conversion
K2 = 1          #equilibrium constant for H2 adsoprtion
K3 = 1          #equilibrium constant for H2O dissociation
p0 = {'CO2': 1, 'H2': 1, 'CO': 0.001, 'H2O': 0.001} #partial pressures
p = p0                                              #for iterative pressure calculations
K = p['CO']*p['H2O']/p['CO2']/p['H2']               #overall equilibrium constant
rate = k1*L0*p0['CO2']/(1 + (K2**0.5)*p0['H2'])     #RWGS rxn rate

#include in iteration:
#rxn rate, in terms of initial partial pressures of CO2 and H2. variable = X
rate = k1*L0*p0['CO2']*(p0['H2']*((1-X)**2) - p0['CO2']*(X**2)/K)
rate = rate/(p0['H2']*(1-X) + (K2**0.5)*p0['H2']*1.5*(((1-X)**1.5)) + p0['CO2']*X/K2/K3)

K = p['CO']*p['H2O']/p['CO2']/p['H2']

print('This is the rate:' + str(rate))
print('This is K:' + str(K))
