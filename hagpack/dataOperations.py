# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 11:11:31 2015

@author: htelg
"""
import numpy as np



def _moving_average(a, n=5) :
    out = np.zeros(a.shape)
    out[:] = np.nan
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    avg = ret[n - 1:] / n
    return avg
    
def moving_average(x,y,
                   windows = [20,100],
                   changeAt = [0.6],
                   ):
    d = x
    amp = y
    
    daOUT = np.array([])
    ampaOUT = np.array([])
    arg_change_last = 0
    for e,ww in enumerate(windows):
        da = _moving_average(d,ww)
        ampa = _moving_average(amp,ww)
        if e > 0:
            arg_change_last = closest(da,changeAt[e-1])
        if e == len(changeAt):
            arg_change = -1
        else:
            arg_change = closest(da,changeAt[e])
#     print arg_change
        da = da[arg_change_last:arg_change]
        ampa = ampa[arg_change_last:arg_change]
        daOUT = np.concatenate((daOUT,da))
        ampaOUT = np.concatenate((ampaOUT, ampa))

    return daOUT,ampaOUT
def closest(a,v):
    return np.abs(a-v).argmin()