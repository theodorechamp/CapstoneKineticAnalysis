

def conversion(t, y0, p):
    numvariables = 4
    y = {'conv':0, 'nvac':0, 'nh2':0, 'nco': 0}
    dadt = p.conversion / p.tspan * p.tstep #maybe change later
    Yco = (p.n * p.tau / p.Vr * dadt
    y['conv'] = p.k1 * Yco**p.gamma1 * (1-y0[1])**p.n1 +\
    p.k3 * Yco**p.gamma3 * (1-y0[1])**p.n1
    y['nvac'] = - y['conv'] * (3 * p.m/p.mw - p.ninitoxy)
    y['nh2'] = - p.x * y['nvac']
    y['nco'] = -(1-p.x) * y['nvac']
