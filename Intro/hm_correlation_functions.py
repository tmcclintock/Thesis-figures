import numpy as np
import cluster_toolkit as ct
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family='serif', size=20)

Om = 0.3
h = 0.7

from colossus.halo import concentration
from colossus.cosmology import cosmology
cos = {'flat':True,'H0':h*100,'Om0':Om,'Ob0':0.05,'sigma8':0.92,'ns':0.96}
cosmology.addCosmology('fiducial', cos)
cosmology.setCosmology('fiducial')

z = 0
r = np.logspace(-2, np.log10(200), 1000)
k = np.loadtxt("datafiles/knl.txt")
P = np.loadtxt("datafiles/pnl.txt")
klin = np.loadtxt("datafiles/klin.txt")
Plin = np.loadtxt("datafiles/plin.txt")
M = 1e14*h #1e14 Msun
c = concentration.concentration(M,'200m',z=z,model='diemer15')

xi_nfw   = ct.xi.xi_nfw_at_R(r, M, c, Om)
xi_mm    = ct.xi.xi_mm_at_R(r, k, P)
bias = ct.bias.bias_at_M(M, klin, Plin, Om)
xi_2halo = ct.xi.xi_2halo(bias, xi_mm)
xi_hm    = ct.xi.xi_hm(xi_nfw, xi_2halo)

fig, ax = plt.subplots(1)

ax.loglog(r/h, xi_nfw, c='r', label=r"NFW")
ax.loglog(r/h, xi_mm, c='b', label=r"Matter")
ax.loglog(r/h, xi_2halo, c='k', ls='--', label=r"$b(M)\xi_{mm}$")

ax.legend(frameon=False)
ax.set_xlabel(r"$r\ [{\rm Mpc}]$")
ax.set_ylabel(r"Correlation function $\xi$")
ax.set_xlim(min(r/h), max(r/h))

#ax.set_yscale('symlog')
#ax.axhline(0, c='k')
#fig.savefig("correlation_function_figure.pdf", bbox_inches='tight')
plt.show()
