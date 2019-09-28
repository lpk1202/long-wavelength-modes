import numpy as np 
from collections import deque
xi = deque()
xj = deque()
xk = deque()

sum_index = 0
for i in np.arange(2,1278):
	print(i)
	for j in np.arange(2,1278):
		for k in np.arange(2,1278):
			if (np.sqrt((2*i-1280)**2+(2*j-1280)**2+(2*k-1280)**2)>38) and (np.sqrt((2*i-1280)**2+(2*j-1280)**2+(2*k-1280)**2)<470):
				xi.append(i)
				xj.append(j)
				xk.append(k)
				sum_index = sum_index+1

print(sum_index)
np.savetxt('data/xi.txt',xi)
np.savetxt('data/xj.txt',xj)
np.savetxt('data/xk.txt',xk)
