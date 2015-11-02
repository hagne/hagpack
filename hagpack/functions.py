import numpy as np

def gauss(x,amp,pos,sigma):
    """ amp,x_center,sigma
        x: array"""
    return 1/(sigma * np.sqrt(2 * np.pi)) * amp * np.exp(-(x-pos)**2/(2.*sigma**2))