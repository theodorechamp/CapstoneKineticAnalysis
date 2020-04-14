import matplotlib.pyplot as plt
import numpy as numpy
import math
from scipy.integrate import solve_ivp

#total mass flow rate (for mass flux in Ergun equation)
def mflow(n, MW):
    mflowTot = 0; #
    for i in range(len(MW)):
        mflowTot = mflowTot + n[i]*MW[i] #[=] g hr**-1
    mflowTot = mflowTot/1000 #[=] kg hr**-1
    return(mflowTot)

#Ergun equation
def ergun(p, MW, T, R, mu, mflowTot, D):
    rho = 0 #total density of the gas mixture
    for i in range(len(p)):
        rho = rho + p[i]*MW[i]/(T*R) #[=] g L**-1 = kg m**-3
    phi = 0.3 #void fraction
    Dp = 0.000385 #[=] m. Particle diameter. Avg of range reported by Gines et al., 1997.
    Ac = math.pi*(D**2)/4 #[=] m**2. cross sectional area
    G0 = mflowTot/Ac/3600 #kg s**-1 m**-2. mass flux
    #differential pressure change. want all units in SI.
    dPdV = -G0*(1-phi)/rho/Ac/Dp/(phi**3)*(150*mu*(1-phi)/Dp + 7/4*G0)
    #conversion from kg m**-4 hr**-2 to atm L**-1
    dPdV = dPdV/(1000**3)/101325
    return(dPdV)

# solve_ivp syntax: f(t, y, c)
# our y = [nCO2, nH2, nCO, nH2O, Ptot]
def reaction(V, y, C):
    # unpacking C
    P0CO2, P0H2, dV, MW, nDict, T, R, mu, mflowTot, D, Vtot = C
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
    Ptot = y[4]

    #kinetic parameters for rate. taken from results of Gines et al., 1997
    #k1L0 = (rate of fwd CO2 dissoc.)*(conc. of empty catalyst sites)
    k1L0 = 0.1287 #units = mol/h/gcat/atm**2
    K2 = 0.039 #equilibrium constant for H2 adsoprtion
    K3 = 0.00000101 #equilibrium constant for H2O dissociation
    K = p[2]*p[3]/p[0]/p[1]
    X = (P0CO2 - p[0])/P0CO2
    #capping conversion at 1.0
    if X >= 1:
        X = 1
    # print('V: ' + str(V))
    # print('K: ' + str(K))
    # print('X: ' + str(X))

    #calculating eta = parameter to fix erroneous rate
    mCat = 3950*Vtot #cat. loading = Al2O3 density multiplied by vol of 1 tube
    W = mCat/nDict['CO2']
    eta = 2464.5*W**2 - 244020*W + 7439600
    #rate [=] mol gCat**-1  hr**-1
    r = k1L0*P0CO2*(P0H2*((1-X)**2))*eta
    r = r/(P0H2*(1-X) + (K2**0.5)*P0H2*1.5*(((1-X)**1.5)) + P0CO2*X/K2/K3)
    # print('rate, V: ' + str(r) + '   ' + str(V))

    # dn/dV = dn/dt/Vdot = dn/dt/ndot/R/T*Ptot = nu*r*gCat*Ptot/(ndot*R*T)
    # first four entries are dn/dV, 5th is dPtot/dV
    dydV = [-1*r*mCat*Ptot/(ntot*R*T), -1*r*mCat*Ptot/(ntot*R*T),
            1*r*mCat*Ptot/(ntot*R*T), 1*r*mCat*Ptot/(ntot*R*T),
            ergun(p, MW, T, R, mu, mflowTot, D)]
    # print('dydV: ' + str(dydV))
    return(dydV)

def Reynolds(p, T, R, n, MW, mu, D):
    # rho = 0 #total density of the gas mixture
    # for i in range(len(p)):
    #     rho = rho + p[i]*MW[i]/(T*R) #g/L = kg m**-3
    mdot = mflow(n, MW) #total mass flow
    mdot = mdot/3600 #conversion to [=] kg s**-1
    Ac = math.pi*(D**2)/4
    Re = mdot*D/mu/Ac #Reynolds number
    return(Re)

def tubes(sol, MWDict):
    #total amt of CO produced in one tube.
    FCO = sol.y[2][-1] #[=] mol CO hr**-1
    BL = 1000000/3 #battery limit = 1 million mt/yr syngas. one third is CO.
    BL = BL/330/24 #conversion to [=] mt/hr
    BL = BL*1000000/MWDict['CO'] #conversion to [=] mol CO hr**-1
    Ntubes = BL/FCO
    return(Ntubes)

def economics():
    #maybe we want to calculate: price of catalyst,
    #capital cost of tubes,
    #cost of energy
    #cost of reactant gases
    #return on recycled CO2 and H2?
    return()

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

def main():
    #constant parameters
    R = 0.08205 #gas constant. [=] L atm K**-1 mol**-1
    T = 523 #isothermal temperature. [=] K
    Hrxn = 42.1 #[=] kJ mol**-1
    Ptot = 44*1.01325 # feed pressure. units converted from bar to atm
    mu = 3.6*0.01849/3600 #3.6 = unit conversion. Value from HYSYS. [=] kg m**-1 s**-1
    #reactor dimensions
    D = 0.0254 #tube diameter. determined from Aboudheir et al. [=] m
    D = 0.05
    # D = 0.1
    L = 1 #tube length. [=] m
    Vtot = 1000*L*math.pi*(D**2)/4 #total volume of 1 tube. units converted to L
    #molecular weight dictionary for all species
    MWDict = {'CO2': 44.008, 'H2': 1.008, 'CO': 28.008, 'H2O': 18.016}
    MW = list(MWDict.values()) #create list instead of dictionary for calculations

    #feed molar flowrates of each species. [=] mol/hr
    ratio = 6
    nDict = {'CO2': 360.7, 'H2': 1, 'CO': 0.0001, 'H2O': 0.0001}
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
            'H2O': Ptot*yDict['H2O']}
    #extract initial CO2 and H2 pressures (needed for rate calculations).
    P0CO2 = pDict['CO2']
    P0H2 = pDict['H2']
    p = list(pDict.values()) #create list instead of dictionary for calculations

    # V will be our t for solve_ivp
    V = numpy.linspace(0, Vtot, 101)
    #list of constant args needed for f(t, y, C). dV = V[1]
    C = [P0CO2, P0H2, V[1], MW, nDict, T, R, mu, mflowTot, D, Vtot]
    # y0 = [nCO20, nH20, nCO0, nH2O0, P0] = molar flowrates, and total pressure.
    y0 = [ nDict['CO2'], nDict['H2'], nDict['CO'], nDict['H2O'], Ptot]
    sol = solve_ivp(lambda t, y: reaction(t,y,C), [V[0], V[-1]], y0, t_eval = V)
    # print('sol.y[0]: ' + str(sol.y[0]))

    #correcting erroneous pressure drop
    dP = [0]*(len(V)-1)
    for i in range(len(V) - 1):
        dP[i] = (sol.y[4][i+1] - sol.y[4][i])*10**6
    #calculating outlet pressure
    Pfinal = Ptot
    for i in range(len(V) - 1):
        Pfinal = Pfinal + dP[i]
    print('Pfinal: ' + str(Pfinal))

    W = 3950*Vtot
    FCO2 = nDict['CO2']
    print('W/FCO2: ' + str(W/FCO2))

    #final conversion
    XCO2 = (nDict['CO2'] - sol.y[0][-1])/nDict['CO2']
    print('XCO2: ' + str(XCO2))

    print('mflowTot: ' + str(mflowTot))
    Re = Reynolds(p, T, R, n, MW, mu, D)
    print('Re: ' + str(Re))

    print('COfinal: ' + str(sol.y[2][-1]))
    Ntubes = tubes(sol, MWDict)
    print('Ntubes: ' + str(Ntubes))

    #total heat
    heat = Hrxn*sol.y[2][-1]*1000/3600 #[=] Watts
    print('Watts: ' + str(heat))

    plotdata(V, sol)

if __name__ == '__main__':
    main()
