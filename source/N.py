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
h=0.6726

k = np.loadtxt('data/k.txt')
pk = np.loadtxt('data/pk.txt')

matter_power_lin = interpolate.interp1d(k, pk)

k_l = 1e-5 # In Mpc^{_1}
print(matter_power_lin(k_l))
def F_2(ks,cos_theta,kl):
	return (5+2*cos_theta*cos_theta)/7-0.5*cos_theta*(ks/kl+kl/ks)

def norm_ks_prime(ks,cos_theta,kl):
	return np.sqrt(ks*ks-2*ks*kl*cos_theta+kl*kl)

def cos_theta_couterpart(ks,cos_theta,kl):
	return (kl-kl*cos_theta)/norm_ks_prime(ks,cos_theta,kl)

def F_2_couterpart(ks,cos_theta,kl):
	return (5+2*cos_theta_couterpart(ks,cos_theta,kl)*cos_theta_couterpart(ks,cos_theta,kl))/7-0.5*cos_theta_couterpart(ks,cos_theta,kl)*(norm_ks_prime(ks,cos_theta,kl)/kl+kl/norm_ks_prime(ks,cos_theta,kl))

def part1(lnks,cos_theta,kl):
	ks = np.exp(lnks)
	return (ks*ks*ks*F_2(ks,cos_theta,kl)*F_2(ks,cos_theta,kl)*matter_power_lin(ks)/(matter_power_lin(norm_ks_prime(ks,cos_theta,kl))))/(4.*np.pi*np.pi)

def part2(lnks,cos_theta,kl):
	ks = np.exp(lnks)
	return (ks*ks*ks*F_2(ks,cos_theta,kl)*F_2_couterpart(ks,cos_theta,kl))/(4.*np.pi*np.pi)

def I(kl): 
	# return 1/(integrate.dblquad(lambda cos_theta, lnks: part1(lnks,cos_theta,kl)+part2(lnks,cos_theta,kl),np.log(kl), np.log(2*kl), lambda x: -1, lambda x: np.exp(x)/(2*kl))[0]\
		return 	1/(integrate.dblquad(lambda cos_theta, lnks: part1(lnks,cos_theta,kl)+part2(lnks,cos_theta,kl),np.log(0.01), np.log(0.15), lambda x: -1, lambda x:1 )[0])


print(I(k_l))
k_plot=k[0:(len(k)-5)]
plt.figure(figsize=(14, 6))
plt.semilogx(k_plot,part1(np.log(k_plot),0.1,k_l), label = 'integrand')
plt.xlabel('$k$')
plt.ylabel('$integrand$')
plt.legend()
plt.savefig('Integrand.png')
plt.show()
