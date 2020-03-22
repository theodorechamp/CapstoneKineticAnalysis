import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import json
from fun_conversionoxidation import *

def execODE(p):
    nvac = (4-p['delta'][0])*p['mH'][0]*p['MWH'][0]

    y0 = {
        'conv': 0,
        'nvac': nvac,
        'nh2': 0,
        'nco': 0
    }
    sol = solve_ivp(conversion, y0, [0, 10],  args = p)
