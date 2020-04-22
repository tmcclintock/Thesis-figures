import numpy as np
from classy import Class
import matplotlib.pyplot as plt
import vapeplot
plt.rc("text", usetex=True)
plt.rc("font", family='serif', size=20)
#cmap = plt.get_cmap("seismic")
cmap = vapeplot.cmap('mallsoft')
colors = [cmap(ci) for ci in np.linspace(1,0,3)]
al = 0.3

G = 4.51715e-48 #Newton's gravitional constant in Mpc^3/s^2/Solar Mass
Mpcperkm = 3.240779289664256e-20 #Mpc/km; used to convert H0 to s^-1
rhocrit = 3.*(Mpcperkm*100)**2/(8*np.pi*G) #Msun h^2/Mpc^3

Omega_b = 0.05 #NOTE: use Omega, not omega, since they mean different things
Omega_m = 0.3
rhom = rhocrit*Omega_m
Omega_cdm = Omega_m - Omega_b
h = 0.7 #This is H0/100
A_s = 2.1e-9 #Power spectrum amplitude
n_s = 0.96 #Power spectrum index
k_max = 0.99
params = {
    'output': 'mPk',
    'non linear':'halofit',
    'h':h,
    #'A_s':A_s,
    'sigma8': 0.9,
    'n_s':n_s,
    #'w0_fld': -.1,
    #'wa_fld': 0.,
    'Omega_b':Omega_b,
    'Omega_cdm':Omega_cdm,
    #'Omega_fld': 1.-Omega_m,
    'P_k_max_1/Mpc':k_max,
    'z_max_pk':10. #Default
    }

k = np.logspace(-4, np.log10(k_max), num=1000) #Mpc^-1
z = 0.

cosmo = Class()
cosmo.set(params)
cosmo.compute()


def calc(params, z, reset=False):
    if reset:
        cosmo.set(params)
        cosmo.compute()
    Pmm  = np.array([cosmo.pk(ki, z) for ki in k]) 
    Plin = np.array([cosmo.pk_lin(ki, z) for ki in k])
    return Plin

if __name__ == "__main__":
    fig, axes = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(12,6))
    for z, i in zip([0., 0.5, 1.], range(3)):
        ax = axes[1]
        Plin = calc(params, z, reset=True)
        ax.loglog(k, Plin, c=colors[i], label=r"$z=%.1f$"%z)

    ax.legend(frameon=False)
    
    #Now the cosmo models
    cmap = vapeplot.cmap('vaporwave')
    colors = [cmap(ci) for ci in np.linspace(0,1,3)]
    z = 0
    ax = axes[0]
    Plin = calc(params, 0, reset=False)
    ax.loglog(k, Plin, c=colors[0], label=r'$\Omega_m = %.2f,\ \sigma_8=%.2f$'%(0.3, cosmo.sigma8()))
    
    Om = 0.35
    params['Omega_cdm'] = Om - 0.05
    cosmo.set(params)
    cosmo.compute()
    Plin = calc(params, 0, reset=False)
    ax.loglog(k, Plin, c=colors[1], label=r'$\Omega_m = %.2f,\ \sigma_8=%.2f$'%(Om, cosmo.sigma8()))

    Om = 0.3
    params['Omega_cdm'] = Om - 0.05
    params['sigma8'] = 0.98
    cosmo.set(params)
    cosmo.compute()
    Plin = calc(params, 0, reset=False)
    ax.loglog(k, Plin, c=colors[2], label=r'$\Omega_m = %.2f,\ \sigma_8=%.2f$'%(Om, cosmo.sigma8()))

    
    ax.legend(frameon=False)
    
    plt.subplots_adjust(wspace=0)
    axes[0].set_xlabel(r"$k\ [{\rm Mpc}^{-1}]$") 
    axes[1].set_xlabel(r"$k\ [{\rm Mpc}^{-1}]$")
    axes[0].set_ylabel(r"$P(k)\  [{\rm Mpc}^{3}]$")
    axes[0].set_xlim(min(k), max(k))
    axes[1].set_xlim(min(k), max(k))
    fig.savefig("powerspectrum_figure.pdf", bbox_inches='tight')
    plt.show()
