# import packages
from xbout import open_boutdataset
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# open dataset
ds = open_boutdataset('BOUT.dmp.*.nc', inputfilepath='BOUT.inp')

# set size for the figures
plt.rcParams["figure.figsize"] = (7, 4)
plt.rcParams.update({"font.size": 10})

# assign the quantities to variables
Ti = ds['Ti'] # ion temperature
Ni = ds['Ni'] # ion density
vi = ds['vi'] # ion velocity
time = ds['t']  # time
phi = ds['phi'] # potential

# plot potential
plt.figure()
plt.plot(time, np.log(np.abs(phi[:, 1, 32, 0])))
plt.title('Evolution of $\phi$ at x=1, y=32, z=0')
plt.ylabel('$\phi$')
plt.xlabel('t [normalised to $\Omega_i$]')
plt.show()

# fit line to find growth rate
slope, intercept = np.polyfit(time, np.log(np.abs(phi[:,1,32,0])),1)

# print the growth rate
print('The growth rate is =', slope)

