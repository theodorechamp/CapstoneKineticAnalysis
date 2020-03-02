import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import json
from fun_conversionoxidation import *

def execODE(params):
    nvac = params.
    y0 = {
        'conv':0,
        'nvac':.25

    }
    sol = solve_ivp(conversion, [0, 10],  args = params)
