

def conversion(t, y0, p):
    numvariables = 4
    y = {'conv':0, 'nvac':0, 'nh2':0, 'nco': 0}
    dadt = p['conversion'][0] / p['tspan'][0] * p['tstep'][0] #maybe change later
    Yco = (p['n'][0] * p['tau'][0] / p['Vr'][0] * dadt) #change to reactor
    y['conv'] = p['k1'][0] * Yco**p['gamma1'][0] * (1-y0['conv'])**p['n1'][0] +\
    p['k3'][0] * Yco**p['gamma3'][0] * (1-y0['conv'])**p['n1'][0]
    y['nvac'] = - y['conv'] * (4 * p['mH'][0]/p['mw'][0] - p['delta'][0]*p['m'][0]/p['mw'][0])
    y['nh2'] = - p['x'][0] * y['nvac']
    y['nco'] = -(1-p['x'][0]) * y['nvac']
