

def conversion(t, y, p):
    dydt = [5]
    dadt = p.conversion / p.tspan * p.tstep #maybe change later
    concentrationCO = p.Vco / p.Ftot *
    dconvdt = params.k1*
    dydt = -.5*y
