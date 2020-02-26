import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d

#Assume installed from github using "git clone --recursive https://github.com/cmbant/CAMB.git"
#This file is then in the docs folders
camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))
sys.path.insert(0,camb_path)
import camb
from camb import model, initialpower, results
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))
pars = camb.set_params(H0=70, ombh2=0.02156, omch2=0.10094, As=2e-9, ns=0.95, tau=0.055)
data= camb.get_background(pars)

print(data.redshift_at_comoving_radial_distance(3072/0.7))
print(data.comoving_radial_distance(1.2,tol=0.0001)*0.7)

x = np.loadtxt('../../../start_now/nbodykit/fastpm-python/Example_halo/data/x_s.txt')
y = np.loadtxt('../../../start_now/nbodykit/fastpm-python/Example_halo/data/y_s.txt')
z = np.loadtxt('../../../start_now/nbodykit/fastpm-python/Example_halo/data/z_s.txt')

r=np.sqrt(x*x+y*y+z*z)
z=np.zeros(len(r))
lll = len(r)
for i in np.arange(lll):
	print(i/lll)
	z[i] = data.redshift_at_comoving_radial_distance(r[i]/0.7)

np.savetxt('z.txt',z)

