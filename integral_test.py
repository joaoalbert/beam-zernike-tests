from scipy.integrate import simps
import numpy as np
import zernike_class as zc


'''
To get the right result in Wolfram Alpha:
integrate 1/sqrt(1-x^2-y^2) dxdy for x=-0.1..0.1, y=-0.1..0.1
'''

# Wolfram answers:
#
# duv=0.01   -   S = 0.000400013
# duv=0.03   -   S = 0.00360108
# duv=0.1    -   S = 0.0401343
# duv=0.4    -   S = 0.678638


duv = 0.1
grid_lim = [-duv,-duv,duv,duv]
Npoints = [101,101]
grid_center = [0,0]
cols = np.sqrt(1/2) + np.zeros((Npoints[0]*Npoints[1],4))

uv_area = (duv*2)**2
print("~uv^2 =", uv_area)



# f(x) unitary: f(x)=1
frequency = 0
beam = zc.Beam_Data(cols, grid_lim, grid_center, Npoints, frequency)
data = beam.co.reshape(Npoints)
beam.generate_grid(fast=True)
thetas_matrix = (beam.thetaphi_grid[:,0]).reshape(Npoints)
#phis_matrix   = (beam.thetaphi_grid[:,1]).reshape(Npoints)
uv = (beam.uv_grid).reshape((Npoints[0],Npoints[1],2))


# f(x) azimuthal gaussian e^(-theta^2)
#data = np.exp(-thetas_matrix**2/(2*0.01**2))/np.sqrt(2*np.pi)/0.01


jacob = 1 / np.cos(thetas_matrix)
#jacob = 1 / np.sqrt(1-uv[:,:,0]**2-uv[:,:,1]**2)
data = data*jacob


I_u = [simps(data[i,:]**2, np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
I_uv =  simps(I_u, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))

print("S_int =",I_uv)



# Monte Carlo

data_max = np.max(data)
# Number of tests per point. Total: N*Npoints[0]*Npoints[1]
N = 10
total = N*Npoints[0]*Npoints[1]
valid = 0

samples = data_max*np.random.rand(Npoints[0],Npoints[1],N)

for i in range(Npoints[0]):
	for j in range(Npoints[1]):
		
		y = data[i][j]
		valid += np.sum(samples[i,j,:]<=y)

print("S_MC  =", valid/total * data_max*uv_area)




