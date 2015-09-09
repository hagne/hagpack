from pint import UnitRegistry
ureg = UnitRegistry()


def mbar2torr(mbar):
    bar = mbar * 1e-3 * ureg.bar
    torr = bar.to(ureg.torr)
    return torr.magnitude
    
def torr2mbar(torr):
    torr = torr * ureg.torr
    mbar = torr.to(ureg.bar) * 1e3
    return mbar.magnitude
    
def ccs2lm(ccs):    
    """ converts cubic centimeters per second to liters per minute"""
    ccs = ccs * ureg.cc / ureg.second
    lm = ccs.to(ureg.liter / ureg.minute)
    return lm.magnitude
    
def lm2ccs(lm):
    """ converts liters per minute to cubic centimeters per second"""
    lm = lm * ureg.liter / ureg.minute
    ccs = lm.to(ureg.cc / ureg.second)
    return ccs.magnitude