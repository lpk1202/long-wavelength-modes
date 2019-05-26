import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
from scipy import interpolate
from sympy.physics.wigner import wigner_3j
import numpy as np
from scipy import interpolate

filename = ''
font = {#'family' : 'normal',
                #'weight' : 'bold',
                'size'  : 18}

plt.rc('font', **font)
#############################################################################################################################################################################################################
### DATA from Cosmosis, DEMO1. I didn't change any of the input values (h=0.6726, omega_m=0.3141, omega_b=0.04, omega_k=0. yhe=0.245341, tau=0.08, n_s=0.96, A_s=2.1e-9, k_s=0.05, w=-1.0.) #################
#############################################################################################################################################################################################################

h=0.6726

k = np.loadtxt('data/k.txt')
pk = np.loadtxt('data/pk.txt')
#############################################################################################################################################################################################################
### pk is the linear power spectrum in Mpc3, k wavevector in Mpc-1; Hubble constant h has been obsorbed. ####################################################################################################
#############################################################################################################################################################################################################

matter_power_lin = interpolate.interp1d(k, pk)

def F_2(ks,cos_theta,kl): # F_2(-ks,ks+ks') in text
	return (5+2*cos_theta*cos_theta)/7-0.5*cos_theta*(ks/kl+kl/ks)

def norm_ks_prime(ks,cos_theta,kl): # |ks'|=|ks-kl|
	return np.sqrt(ks*ks-2*ks*kl*cos_theta+kl*kl)

def cos_theta_couterpart(ks,cos_theta,kl):
	return (kl-ks*cos_theta)/norm_ks_prime(ks,cos_theta,kl)

def F_2_couterpart(ks,cos_theta,kl): # F_2(-ks',ks+ks') in text
	return (5+2*cos_theta_couterpart(ks,cos_theta,kl)*cos_theta_couterpart(ks,cos_theta,kl))/7-0.5*cos_theta_couterpart(ks,cos_theta,kl)*(norm_ks_prime(ks,cos_theta,kl)/kl+kl/norm_ks_prime(ks,cos_theta,kl))

def part1(lnks,cos_theta,kl):
	ks = np.exp(lnks)
	return (ks*ks*F_2(ks,cos_theta,kl)*F_2(ks,cos_theta,kl)*matter_power_lin(ks)/(2*matter_power_lin(norm_ks_prime(ks,cos_theta,kl))))/(4.*np.pi*np.pi)

def part2(lnks,cos_theta,kl):
	ks = np.exp(lnks)
	return (ks*ks*F_2(ks,cos_theta,kl)*F_2_couterpart(ks,cos_theta,kl))/(4.*np.pi*np.pi)

def part3(lnks,cos_theta,kl):
	ks = np.exp(lnks)
	return (ks*ks*F_2_couterpart(ks,cos_theta,kl)*F_2_couterpart(ks,cos_theta,kl)*matter_power_lin(norm_ks_prime(ks,cos_theta,kl))/(2*matter_power_lin(ks)))/(4.*np.pi*np.pi)

def I(kl): 
	# return 1/(integrate.dblquad(lambda cos_theta, lnks: part3(lnks,cos_theta,kl), np.log(kl), np.log(2*kl), lambda x: -1, lambda x: np.exp(x)/(2*kl))[0]\
	return 1/(integrate.dblquad(lambda cos_theta, lnks: part1(lnks,cos_theta,kl)+part2(lnks,cos_theta,kl)+part3(lnks,cos_theta,kl),np.log(0.01), np.log(0.15), lambda x: -1, lambda x:1 )[0])

Signal = np.zeros(100)
Noise = np.zeros(100)
X = np.logspace(-4,np.log10(3e-3),100)

for i in np.arange(len(X)):
	Signal[i] = matter_power_lin(X[i])
	Noise[i] = I(X[i])
	print(i,' ',X[i],' ',Signal[i],' ', Noise[i])


# k_plot=k[0:(len(k)-5)]
plt.figure(figsize=(14, 6))
plt.loglog(X, Signal, label = '1')
plt.loglog(X, Noise, label = '1')
plt.xlabel(r'$k$ in $\rm Mpc^{-1}$')
plt.ylabel(r'Power Spectrum and Its Variance in $\rm Mpc^{-3}$')
plt.legend([r'$P(k)$',r'$N(k)$'])
plt.savefig('StoN.png')
plt.show()
