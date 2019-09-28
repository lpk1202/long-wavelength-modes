import sys
sys.path.insert(0, '../')

from nbodykit.lab import *
import matplotlib.pyplot as plt
from nbodykit import setup_logging
from fastpm.nbkit import FastPMCatalogSource
from pmesh.pm import ParticleMesh, RealField, ComplexField
from fastpm import core
import numpy as np 
from fastpm.core import leapfrog
from nbodykit.lab import BigFileMesh
from pmesh.pm import ParticleMesh
from nbodykit.cosmology import Planck15
from kdcount import KDTree
from scipy import interpolate

from nbodykit.source.catalog import ArrayCatalog

# generate the fake data

x = np.loadtxt('../../../nb/data/x.txt')
y = np.loadtxt('../../../nb/data/y.txt')
z = np.loadtxt('../../../nb/data/z.txt')
mass = np.loadtxt('../../../nb/data/m_c200.txt')
length = len(x)
data = np.empty(length, dtype=[('Position', ('f8', 3)), ('Mass', 'f8')])
pos = np.zeros((length,3))

for i in np.arange(length):
	pos[i][0] = x[i]
	pos[i][1] = y[i]
	pos[i][2] = z[i]
data['Position'] = pos
data['Mass'] = mass

# initialize the catalog
f = ArrayCatalog(data)

print(f)
print("columns = ", f.columns) # default Weight,Selection also present
print("total size = ", f.csize)


f = ArrayCatalog({'Position' : data['Position'], 'Mass' : data['Mass'] })

print(f)
print("columns = ", f.columns) # default Weight,Selection also present
print("total size = ", f.csize)

# convert to a MeshSource, using CIC interpolation on 1280^3 mesh
mesh = f.to_mesh(window='cic', Nmesh=1280, BoxSize = 4000.0, compensated=True,interlaced=True)

rfield = mesh.compute()
cfield = (rfield-1).r2c()

cfield_real = np.zeros((1280,1280,641),dtype=float)
cfield_imag = np.zeros((1280,1280,641),dtype=float)
for i in np.arange(1280):
	print(i)
	for j in np.arange(1280):
		for k in np.arange(641):
			cfield_real[i][j][k] = cfield[i][j][k].real
			cfield_imag[i][j][k] = cfield[i][j][k].imag

np.savetxt('data/cfield_real.txt',cfield_real.flatten())
np.savetxt('data/cfield_imag.txt',cfield_imag.flatten())

# cfield2 = np.zeros((1280,1280,1280),dtype=complex)




# for i in np.arange(1280):
# 	for j in np.arange(1280):
# 		for k in np.arange(641):
# 			cfield2[i][j][k] = cfield[i][j][k]

# for i in np.arange(1,1280):
# 	for j in np.arange(1,1280):
# 		for k in np.arange(641,1280):
# 			cfield2[i][j][k] = np.conj(cfield[1280-i][1280-j][1280-k])

# # compute the power, specifying desired linear k-binning
# r = FFTPower(mesh, mode='1d', dk=0.005, kmin=0.005)

# # the result is stored at "power" attribute
# Pk = r.power
# print(Pk)

# # # print the shot noise subtracted P(k)
# # plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'])

# # # format the axes
# # plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
# # plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
# # plt.xlim(0.01, 0.6)
# # plt.savefig('1.png')
# # plt.show()

# h=0.6777
# k_camb = np.loadtxt('data/k_h.txt')
# pk_camb = np.loadtxt('data/p_h.txt')
# matter_power_lin = interpolate.interp1d(k_camb, pk_camb)

# delta_k = np.pi/4000

# res = np.zeros((11,11,11),dtype=complex)
# res_kl = np.zeros((11,11,11),dtype=complex)

# def g(l,m,n,Nlx,Nly,Nlz):
# 	M_kl = np.sqrt(Nlx**2+Nly**2+Nlz**2)
# 	M_ks = np.sqrt((2*l-1280)**2+(2*m-1280)**2+(2*n-1280)**2)
# 	dot_ls = Nlx*(2*l-1280)+Nly*(2*m-1280)+Nlz*(2*n-1280)
# 	M_ks_counterpart = np.sqrt((2*l-1280-Nlx)**2+(2*m-1280-Nly)**2+(2*n-1280-Nlz)**2)
# 	dot_ls_counterpart = -Nlx*(2*l-1280-Nlx)-Nly*(2*m-1280-Nly)-Nlz*(2*n-1280-Nlz)
# 	F2 = 5/7-0.5*dot_ls*(M_kl/M_ks+M_ks/M_kl)/(M_kl*M_ks)+2*(dot_ls)**2/(7*M_kl*M_kl*M_ks*M_ks)
# 	F2_counterpart = 5/7-0.5*dot_ls_counterpart*(M_kl/M_ks_counterpart+M_ks_counterpart/M_kl)/(M_kl*M_ks_counterpart)+2*(dot_ls_counterpart)**2/(7*M_kl*M_kl*M_ks_counterpart*M_ks_counterpart)
# 	return (F2*matter_power_lin(M_ks*delta_k)+F2_counterpart*matter_power_lin(M_ks_counterpart*delta_k))/(2*matter_power_lin(M_ks*delta_k)*matter_power_lin(M_ks_counterpart*delta_k))

# xi = np.loadtxt('data/xi.txt')
# xj = np.loadtxt('data/xj.txt')
# xk = np.loadtxt('data/xk.txt')
# length_x = len(xi)

# for N_kl_x in np.arange(-10,12,2): #have to be even numbers 
# 	for N_kl_y in np.arange(-10,12,2):
# 		for N_kl_z in np.arange(-10,12,2):
# 			Int = 0.
# 			for i in np.arange(length_x):
# 				print(i,' ',N_kl_x,' ',N_kl_y,' ',N_kl_z)
# 				Int  = Int + g(xi[i],xj[i],xk[i],N_kl_x,N_kl_y,N_kl_z)*cfield2[xi[i]][xj[i]][xk[i]]*cfield2[1280+int(N_kl_x/2)-xi[i]][1280+int(N_kl_y/2)-xj[i]][1280+int(N_kl_z/2)-xk[i]]
			 
# 			res[int(N_kl_x/2)+5][int(N_kl_y/2)+5][int(N_kl_z/2)+5] = Int
# 			res_kl[int(N_kl_x/2)+5][int(N_kl_y/2)+5][int(N_kl_z/2)+5] = cfield2[640+int(N_kl_x/2)][640+int(N_kl_y/2)][640+int(N_kl_z/2)]
# 			print(res[int(N_kl_x/2)+5][int(N_kl_y/2)+5][int(N_kl_z/2)+5]*4*np.pi*((upperlimit/h)**3-(lowerlimit/h)**3)/(3*sum_index*(2*np.pi)**3),'akdnasoifio',res_kl[int(N_kl_x/2)+5][int(N_kl_y/2)+5][int(N_kl_z/2)+5])

# sum_index = 54330094

# upperlimit = 470*np.pi/(4000/h)
# lowerlimit = 38*np.pi/(4000/h)
# res = res*4*np.pi*((upperlimit/h)**3-(lowerlimit/h)**3)/(3*sum_index*(2*np.pi)**3)
# np.savetxt('data/res.txt',res)
# np.savetxt('data/res_kl.txt',res_kl)



