import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
import vapeplot
plt.rc("text", usetex=True)
plt.rc("font", family='serif', size=20)
#cmap = plt.get_cmap("seismic")
cmap = vapeplot.cmap('mallsoft')
colors = [cmap(ci) for ci in np.linspace(1,0,3)]
al = 0.3


def calc(zs, Om, w=None):
    if w == None:
        p6=ccl.Parameters(Omega_c=Om-0.05, Omega_b=0.05, h=0.7, sigma8=0.9, n_s=0.96)
    else:
        p6=ccl.Parameters(Omega_c=Om-0.05, Omega_b=0.05, h=0.7, sigma8=0.9, n_s=0.96, w0=w)
    cosmo = ccl.Cosmology(p6)

    sfs = 1./(1+zs)
    return ccl.background.growth_factor(cosmo, sfs)

if __name__ == "__main__":
    fig, ax = plt.subplots(1)
    sfs = np.logspace(-4, 0)
    zs = 1./sfs - 1
    g = calc(zs, 0.3)
    ax.plot(sfs, g, label=r"$\Omega_m = %.1f\ w=%.0f$"%(0.3, -1))

    g = calc(zs, 0.3, -0.7)
    ax.plot(sfs, g, label=r"$\Omega_m = %.1f\ w=%.1f$"%(0.3, -0.7))

    g = calc(zs, 0.25)
    ax.plot(sfs, g, label=r"$\Omega_m = %.1f\ w=%.0f$"%(0.25, -1))

    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlabel(r"Scale factor $a$")
    ax.set_ylabel(r"Growth function $g(a)/g(1)$")


    ax.plot(sfs, sfs, label=r"$Eds, \Omega_m=1$")
    ax.legend(frameon=False, fontsize=12)
    fig.savefig("growth_figure.pdf", bbox_inches='tight')
    plt.show()
