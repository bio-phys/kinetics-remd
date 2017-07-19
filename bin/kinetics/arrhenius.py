__author__ = 'lustelzl'

from ala_kinetics import *
from scipy.optimize import curve_fit

def arrhenius(T, Ea, A):
    R=8.3144621
    return -Ea/R * (1.0/T) + np.log(A)
