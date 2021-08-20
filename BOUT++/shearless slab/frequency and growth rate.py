# import packages
from xbout import open_boutdataset
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Model

# open dataset
ds = open_boutdataset('BOUT.dmp.*.nc', inputfilepath='BOUT.inp')

# set size for the figures
plt.rcParams["figure.figsize"] = (7, 4)
plt.rcParams.update({"font.size": 10})

# assign variables
Ti = ds['Ti'] # assign ion temperature
Ni = ds['Ni'] # assign ion density
time = ds['t'] # assign time
phi = ds['phi'] # assign phi


# isolate region for gamma
time_region_gamma = time[300:]
phi_region = phi[300:, 1, 16, 4]

# isolate region for omega
Ni_region = Ni[500:, 1, 16, 4]
time_region_omega = time[500:]

# guesses for (cosine) fit
amplit = 0.000003 # amplitude
freq = 0.8 # frequency
c = 1 # phase

def linear_fit(x, a, b):
    return a*x+b

# linear fit for the growth rate
popt, pcov = curve_fit(linear_fit, time_region_gamma,np.log(np.abs(phi_region)))

slope = popt[0]
gamma_st_dev = np.sqrt(pcov[0][0])

# plot phi against time and the fit
plt.plot(time_region_gamma, np.log(np.abs(phi_region)), label='data')
plt.plot(time_region_gamma, linear_fit(time_region_gamma, *popt),'--', label='fitting')
plt.xlabel('t')
plt.ylabel('phi')
#plt.yscale('log')
plt.legend()
plt.show()

print('growth rate  = ', popt[0], '+/-', np.sqrt(pcov[0][0]))

gamma = slope
# define a cosine function for fitting
def cos_func(x, amplitude, frequency, phase):

    return amplitude*np.cos(frequency*x + phase)
    
# fit with lmfit
gmodel = Model(cos_func)
result = gmodel.fit(Ni_region/np.exp(slope*time_region_omega), x=time_region_omega, amplitude=amplit, frequency=freq, phase=c)
print(result.fit_report())

plt.plot(time_region_omega, Ni_region/np.exp(slope*time_region_omega), label='data')
plt.plot(time_region_omega, result.best_fit,'--', label='fitting')
plt.ylabel('Ni/exp(gamma*t)')
plt.xlabel('t')
plt.title('Fitting with lmfit')
plt.legend()
plt.show()
