import numpy as np

def ReynoldsNoInPipe_fromFlowVelocity(v,d,roh = 1.2,mu= 18.27e-6):
    """http://en.wikipedia.org/wiki/Reynolds_number
    v:    \t mean velocity of the fluid (SI units: m/s)
    d:    \t diameter in m
    roh:  \t density of the fluid (kg/m^3) [roh_air(sealevel,T_room) = 1.2]
    mu:   \t dynamic viscosity of the fluid (Pa*s = N*s/m^2 = kg/(m*s)) [mu_air(T=291.15K) = 18.27e-6], mu_air,T=291.15K) = 13.29e-6"""

    ry = (v*d*roh)/mu
    if ry > 2000:
        print 'warning: reynolds number is larger than 2000'
    return ry

def ReynoldsNoInPipe_fromVFlowRate(fl,d,roh = 1.2,mu= 18.27e-6):
    """http://en.wikipedia.org/wiki/Reynolds_number
    fl:   \t volumetric flow rate (cc/s).
    d:    \t diameter in m
    roh:  \t density of the fluid (kg/m^3) [roh_air(sealevel,T_room) = 1.2]
    mu:   \t dynamic viscosity of the fluid (Pa*s = N*s/m^2 = kg/(m*s)) [mu_air(T=291.15K) = 18.27e-6], mu_air,T=291.15K) = 13.29e-6"""
    A = np.pi*(d/2.)**2
    ry = (fl*1e-6*d*roh)/(mu*A)
    if ry > 2000:
        print 'warning: reynolds number is larger than 2000'
    return ry