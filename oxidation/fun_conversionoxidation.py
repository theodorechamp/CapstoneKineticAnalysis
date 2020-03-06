

def conversion(t, y0, p):
    numvariables = 4
    y = {'conv':0, 'nvac':, 'nh2':0, 'nco': 0}
    dadt = p.conversion / p.tspan * p.tstep #maybe change later
    Yco = (p.n * p.tau / p.Vr * dadt) #change to reactor
    y['conv'] = p.k1 * Yco**p.gamma1 * (1-y0['conv'])**p.n1 +\
    p.k3 * Yco**p.gamma3 * (1-y0['conv'])**p.n1
    y['nvac'] = - y['conv'] * (4 * p.mH/p.mw - p.delta*p.m/p.mw)
    y['nh2'] = - p.x * y['nvac']
    y['nco'] = -(1-p.x) * y['nvac']
