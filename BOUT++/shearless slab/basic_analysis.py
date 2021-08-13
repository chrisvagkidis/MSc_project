# import packages
from xbout import open_boutdataset
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# set size for the figures
plt.rcParams["figure.figsize"] = (7, 4)
plt.rcParams.update({"font.size": 10})

# open dataset
ds = open_boutdataset('BOUT.dmp.*.nc', inputfilepath='BOUT.inp')

# plot a specific x point at a specific time for the ion temperature
plt.figure()
ds['Ti'].isel(t=145, x=1).plot()
plt.show()

# plot a specific x point at a specific time for the ion density
plt.figure()
ds['Ni'].isel(t=145, x=1).plot()
plt.show()

# assign the quantities to variables
Ti = ds['Ti'] # ion temperature
Ni = ds['Ni'] # ion density
vi = ds['vi'] # ion velocity
time = ds['t'] # time
phi = ds['phi'] # potential


print('Ti shape = ', Ti.shape)
print('Ni shape = ', Ni.shape)
print('vi shape = ', vi.shape)
print('time shape = ', time.shape)

# create subplots with T_i, N_i, v_i, phi

fig, axs = plt.subplots(2,2, figsize=(15, 7))

# ion temperature
axs[0, 0].plot(time, np.log(np.abs(Ti[:, 1, 32, 0])))
axs[0, 0].set_title('Evolution of $T_i$ at x = 1, y = 32, z = 0')
axs[0, 0].set(xlabel='t', ylabel='$T_i$')

# ion density
axs[0, 1].plot(time,np.log(np.abs(Ni[:, 1, 32, 0])))
axs[0, 1].set_title('Evolution of $N_i$ at x = 1, y = 32, z = 0')
axs[0, 1].set(xlabel='t', ylabel='$N_i$')
#axs[0, 1].semilogy()

# ion velocity
axs[1, 0].plot(time, np.log(np.abs(vi[:, 1, 32, 0])))
axs[1, 0].set_title('Evolution of $v_i$ at x = 1, y = 32, z = 0')
axs[1, 0].set(xlabel='t', ylabel='$v_i$')
#axs[1, 0].semilogy()

# phi
axs[1, 1].plot(time, np.log(np.abs(phi[:, 1, 32, 0])))
axs[1, 1].set_title('Evolution of $\phi$ at x = 1, y = 32, z = 0')
axs[1, 1].set(xlabel='t', ylabel='$\phi$')
#axs[1, 1].set_yscale('log')

plt.show()
