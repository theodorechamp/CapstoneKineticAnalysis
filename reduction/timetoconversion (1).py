
# coding: utf-8

# In[ ]:


import math as m

def timetoconversion(p):
    fint = 0.5*(((1-p.conversion)^(-2))-1) # the integral of the f(alpha) fxn; pressure dependence not needed
    rhs = p.A*exp(-p.Ea/(p.R*p.T)) # right hand side of eqn
    time = fint/rhs # solving for time

