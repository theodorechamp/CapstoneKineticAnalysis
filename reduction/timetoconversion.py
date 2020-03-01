import math as m

def timetoconversion(p):
    fint = -2/(1-p.conversion)^2 + 2 # the integral of the f(alpha) fxn
    pdep = 1 # add in after we find the pressure dependence; may need to take integral of dependence
    rhs = p.A*m.exp(-p.Ea/(p.R*p.T))
    time = (fint/pdep)/rhs
