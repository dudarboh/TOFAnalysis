import numpy as np
import matplotlib.pyplot as plt
c = 299.792458
p = 1.7
l0 = 2050.
t0 = 7.12057
m0 = 493.677

t_lim = l0/c
t = np.arange(t_lim, 10., 0.00001)

l_lim = t0*c
l = np.arange(2000., l_lim, 0.1)

# vs t
# beta = l0/(t*c)
# m = p/beta * np.sqrt(1. - beta**2) * 1000.
# plt.plot(t*1000 - t0*1000, m - m0)
# plt.xlabel(r'$t - t_{true}$, [ps]')

# vs l
beta = l/(t0*c)
m = p/beta * np.sqrt(1. - beta**2) * 1000.
plt.plot(l - l0, m - m0)
plt.xlabel(r'$l - l_{true}$, [mm]')


plt.ylabel(r"$m - m_{true}$, [MeV]")
plt.grid(True)
plt.show()
