import matplotlib.pyplot as plt
import numpy as numpy
import math
from scipy.integrate import solve_ivp

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
        rho = rho + p[i]*MW[i]/(T*R) #g/L = kg/m**3
    phi = 0.3 #void fraction
    Dp = D/8 #[] = m. Particle diameter. Using heuristic Weimer gave in class on 10/10/19 that particle diameter should be 1/8 of tube diameter
    Ac = math.pi*(D**2)/4 #[] = m^2. cross sectional area
    G0 = mflowTot/Ac/3600 #g/s/m^2. mass flux
    mu = 3.6*0.01849/3600 #3.6 = unit conversion. Value from HYSYS. [] = kg m**-1 s**-1
    #differential pressure change. want all units in SI.
    dPdV = -G0*(1-phi)/rho/Ac/Dp/(phi**3)*(150*mu*(1-phi)/Dp + 7/4*G0)
    #convert units from kg m**-4 hr**-2 to atm L**-1
    dPdV = dPdV/(1000**3)/(101325)
    return(dPdV)

# solve_ivp syntax: f(t, y, c)
# our y = [nCO2, nH2, nCO, nH2O, P]
def reaction(V, y, C):
    # unpacking C
    P0CO2, P0H2, dV, MW, T, R, mflowTot, D, Vtot = C
    #cap each species at zero
    for i in range(0,4):
        if y[i] <= 0:
            y[i] = 0
    #calculate p's from y = [nCO2, nH2, nCO, nH2O, Ptot]
    ntot = 0
    for i in range(0,4):
        ntot = ntot + y[i]
    p = [0]*4
    for i in range(0,4):
        p[i] = y[i]/ntot*y[4]
        if y[i] <= 0:
            p[i], y[i] = 0, 0
    #parameters needed to find rate
    #kinetic parameters. taken from results of Gines et al., 1997
    #k1L0 = (rate of fwd CO2 dissoc.)*(conc. of empty catalyst sites)
    k1L0 = 0.1287 #units = mol/h/gcat/atm**2
    K2 = 0.039 #equilibrium constant for H2 adsoprtion
    K3 = 0.00000101 #equilibrium constant for H2O dissociation
    K = p[2]*p[3]/p[0]/p[1]
    X = (P0CO2 - p[0])/P0CO2
    if X >= 1:
        X = 1
    #rate [=] mol gCat**-1  hr**-1
    print('V: ' + str(V))
    print('K: ' + str(K))
    print('X: ' + str(X))
    print('Term1: ' + str(P0H2*((1-X)**2)))
    print('Term2: ' + str(P0CO2*(X**2)/K))
    # print('Term3: ' + str((P0H2*(1-X))))
    # print('Term4: ' + str((K2**0.5)*P0H2*1.5*(((1-X)**1.5))))
    # print('Term5: ' + str(P0CO2*X/K2/K3))
    F = 0.00001 #parameter to fix erroneous rate
    r = k1L0*P0CO2*(P0H2*((1-X)**2)/F - F*P0CO2*(X**2)/K)
    r = r/(P0H2*(1-X) + (K2**0.5)*P0H2*1.5*(((1-X)**1.5)) + P0CO2*X/K2/K3)
    # r = 10
    print('rate: ' + str(r))
    #parameters for calculating dydV
    mCat = 3950*Vtot #cat. loading = Al2O3 density multiplied by vol of 1 tube
    Ptot = y[4]
    # dn/dV = dn/dt/Vdot = dn/dt/ndot/R/T*Ptot = nu*r*gCat*Ptot/(ndot*R*T)
    # first four entries are dn/dV, 5th is dPtot/dV
    # dydV = [-1*r*mCat*dV*Ptot/(ntot*R*T), -1*r*mCat*dV*Ptot/(ntot*R*T),
    #         1*r*mCat*dV*Ptot/(ntot*R*T), 1*r*mCat*dV*Ptot/(ntot*R*T),
    #         ergun(p, MW, T, R, mflowTot, D)]
    dydV = [-1*r*mCat*Ptot/(ntot*R*T), -1*r*mCat*Ptot/(ntot*R*T),
            1*r*mCat*Ptot/(ntot*R*T), 1*r*mCat*Ptot/(ntot*R*T),
            ergun(p, MW, T, R, mflowTot, D)]
    print('dydV: ' + str(dydV))
    return(dydV)

def plotdata(V, sol):
    y = sol.y
    CO2 = y[0]
    H2 = y[1]
    CO = y[2]
    H2O = y[3]
    plt.figure(1)
    plt.plot(V, CO2, 'r', label = 'CO2')
    plt.plot( V, H2, 'g', label = 'H2')
    plt.plot( V, CO, 'm', label = 'CO')
    plt.plot( V, H2O, 'b', label = 'H2O')
    plt.xlabel('Volume (L)')
    plt.ylabel('Flowrate (mol/hr)')
    plt.legend()
    plt.savefig('results.png')
    #plt.show()

def main():
    #constant parameters
    R = 0.08205 #gas constant. [=] L atm K**-1 mol**-1
    T = 523 #isothermal temperature. [=] K
    Ptot = 44*1.01325 # feed pressure. units converted from bar to atm
    #reactor dimensions
    D = 0.0254 #tube diameter. determined from Aboudheir et al. [=] m
    L = 1 #tube length. [=] m
    Vtot = 1000*L*math.pi*(D**2)/4 #total volume of 1 tube. units converted to L

    #molecular weight dictionary for all species
    MWDict = {'CO2': 44.008, 'H2': 1.008, 'CO': 28.008, 'H2O': 18.016}
    MW = list(MWDict.values()) #create list instead of dictionary for calculations

    #feed molar flowrates of each species. [=] mol/hr
    ratio = 3
    nDict = {'CO2': 25, 'H2': 1, 'CO': 0.001, 'H2O': 0.001}
    nDict['H2'] = ratio*nDict['CO2'] - nDict['CO'] - nDict['H2O']
    n = list(nDict.values())
    ntot = 0 #total molar flowrate
    for i in range(len(n)):
        ntot = ntot + n[i]
    mflowTot = mflow(n, MW) #total mass flow
    #feed mole fraction dictionary
    yDict = {'CO2': nDict['CO2']/ntot, 'H2': nDict['H2']/ntot, 'CO': nDict['CO']/ntot, 'H2O': nDict['CO2']/ntot}
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

    #list of constant args needed for f(t, y, C)
    V = numpy.linspace(0, Vtot, 101)
    #dV = V[1]
    C = [P0CO2, P0H2, V[1], MW, T, R, mflowTot, D, Vtot]
    # y0 = [nCO20, nH20, nCO0, nH2O0, P0] = molar flowrates, and total pressure.
    y0 = [ nDict['CO2'], nDict['H2'], nDict['CO'], nDict['H2O'], Ptot]
    #(V, y, P0CO2, P0H2, k1L0, K2, K3, dV, MW, T, R, mflowTot, D, Vtot)
    sol = solve_ivp(lambda t, y: reaction(t,y,C), [V[0], V[-1]], y0, t_eval = V)
    # print('sol: ' + str(sol))
    print('sol.y[0]: ' + str(sol.y[0]))

    #finding total pressure at reactor exit
    dPdVtot = 0
    for i in range(len(V)):
        dPdVtot = dPdVtot + V[1]*sol.y[4][i]
    print('dPdVtot: ' + str(dPdVtot))
    Pfinal = 44*1.01325 - dPdVtot*Vtot
    print('Pfinal: ' + str(Pfinal))

    #final conversion
    XCO2 = (nDict['CO2'] - sol.y[0][-1])/nDict['CO2']
    print('XCO2: ' + str(XCO2))

    plotdata(V, sol)

if __name__ == '__main__':
    main()
