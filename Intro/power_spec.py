import numpy as np
from classy import Class

Omega_b = 0.05 #NOTE: use Omega, not omega, since they mean different things
Omega_m = 0.3
rhom = rhocrit*Omega_m
Omega_cdm = Omega_m - Omega_b
h = 0.7 #This is H0/100
A_s = 2.1e-9 #Power spectrum amplitude
n_s = 0.96 #Power spectrum index
k_max = 1
params = {
    'output': 'mPk',
    'non linear':'halofit',
    'h':h,
    #'A_s':A_s,
    'sigma8': 0.77,
    'n_s':n_s,
    'w0_fld': -.7,
    'wa_fld': 0.,
    'Omega_b':Omega_b,
    'Omega_cdm':Omega_cdm,
    'Omega_fld': 1.-Omega_m,
    'P_k_max_1/Mpc':k_max,
    'z_max_pk':10. #Default
    }

k = np.logspace(-5, np.log10(k_max), base=10, num=1000) #Mpc^-1
z = 0.

def calc(params, z):
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    #Pmm  = np.array([cosmo.pk(ki, z) for ki in k]) 
    Plin = np.array([cosmo.pk_lin(ki, z) for ki in k])
    return Plin

if __name__ == "__main__":
    import matplotlib.pyplot as plt
