import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

def make_histogram(xx_in, dmag, cut_off, min_val = None, cumulative=True):
    xx = xx_in[np.where(xx_in<=cut_off+dmag)]
    #print xx.min(),xx.max()
    if min_val is None:
        min_val=xx.min()-dmag
    i_xx = np.round((xx-min_val)/dmag).astype(int)
    unique_ixx, ct = np.unique(i_xx, return_counts=True)

    if cumulative:
        return unique_ixx*dmag+min_val, ct.astype(float)/float(len(xx_in))
    else:
        return unique_ixx*dmag+min_val, ct.astype(int)

if __name__ == "__main__":

    dtype = np.dtype([('z', float), ('tau', float)])
    data = np.genfromtxt('agn_tau_distribution.txt', dtype=dtype)
    plt.figsize=(30,30)

    tau_renorm = np.log10(data['tau']/(1.0+data['z']))
    tau  = np.log10(data['tau'])
    tau_min = tau.min()
    tau_max = tau.max()
    tau_renorm_min = tau_renorm.min()
    tau_renorm_max = tau_renorm.max()
    tau_min = min(tau_min, tau_renorm_min)
    tau_max = max(tau_max, tau_renorm_max)
    dtau = 0.1
    tau_grid, tau_hist = make_histogram(tau, dtau, tau_max+dtau,
                                        cumulative=False)

    (tau_renorm_grid,
     tau_renorm_hist) = make_histogram(tau_renorm, dtau, tau_max+dtau,
                                       cumulative=False)

    t_l, = plt.plot(tau_grid, tau_hist)
    t_r_l, = plt.plot(tau_renorm_grid, tau_renorm_hist)
    plt.legend([t_l, t_r_l],['$\\tau$', '$\\tau/(1+z)$'], loc=0)
    plt.xlim(0, 5)
    plt.xlabel('$\log(\\tau)$')
    plt.savefig('agn_tau_dist_fig.png')
