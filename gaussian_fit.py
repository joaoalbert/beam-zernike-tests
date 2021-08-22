import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import zernike_fit as zk


#filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/npoints_test/980_1383pts_04uv/spherical_grid.grd"
#filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/npoints_test/980_922pts_005uv/spherical_grid.grd"

#filename = ""
filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/bingo_standard/Retangular_MINUS1110_110_linear_x_980.0/spherical_grid.grd"
#filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/bingo_standard/Retangular_30_MINUS304_linear_x_980.0/spherical_grid.grd"

# *************** CHANGE THE TITLE ACCORDING TO FILE ****************
title = "Horn (-1110,110) - 980MHz "

radius = 0.01

plot_2d = True # if False, will produce surface plots






# *************************************************



def gauss_2d_wrapper(radius=-1):

	def gauss_2d(coords, A, mu_x, mu_y, sigma):
	
		x,y = coords
		x_0, y_0 = x-mu_x, y-mu_y
		
		r2 = x_0**2 + y_0**2
		r02 = radius**2
		
		if radius==-1: return A*np.exp( -r2 /sigma**2 /2 )
		else: return A*np.exp(-r2 /sigma**2 /2) * (r2<=r02) #np.piecewise(r2, [r2<=r02, r2>r02], [A*np.exp(-r2 /sigma**2 /2), 0])
		
	return gauss_2d
	
	
	
def gauss_elliptical_wrapper(radius=-1):

	def gauss_elliptical(coords, A, mu_x, mu_y, sigma_x, sigma_y, alpha):
	
		a = (np.cos(alpha)/sigma_x)**2/2  + (np.sin(alpha)/sigma_y)**2/2
		b = -np.sin(2*alpha)/sigma_x**2/4 + np.sin(2*alpha)/sigma_y**2/4
		c = (np.sin(alpha)/sigma_x)**2/2  + (np.cos(alpha)/sigma_y)**2/2
		
		x,y = coords
		x_0, y_0 = x-mu_x, y-mu_y
		
		r2 = x_0**2 + y_0**2
		r02 = radius**2
		
		if radius==-1: return A*np.exp( a*x_0**2 + 2*b*x_0*y_0 + c*y_0**2 )
		else: return A*np.exp(-( a*x_0**2 + 2*b*x_0*y_0 + c*y_0**2 )) * (r2<=r02)
		
	return gauss_elliptical
	



print("Collecting data...")

cols, grid_lims, grid_centers, Nps, freqs = zk.grasp_spher_grd_file(filename, shift_center=True, verbose=False)	
grid_lim = grid_lims[0]
grid_center = grid_centers[0]
Npoints = Nps[0]
cols = cols[0]

data_0 = np.sqrt(cols[:,0]**2 + cols[:,1]**2)
A, mu_x, mu_y, sigma   = np.max(data_0), 0,0, 0.005 
# grid_center = zk.uv_to_polar(grid_center)
#mu_x, mu_y = center[0]*np.cos(center[1]), center[0]*np.sin(center[1])



print("Making grid...")

uv_grid = [ [u,v] for u in np.linspace(grid_lim[0],grid_lim[2],Npoints[0]) for v in np.linspace(grid_lim[1],grid_lim[3],Npoints[1]) ]
thetaphi_coordinates = zk.polar_uv_grid(uv_grid)

max_point = np.argmax(data_0, axis=None)
grid_center = thetaphi_coordinates[max_point]

coordinates = zk.rotate_coords(thetaphi_coordinates, grid_center, verbose=True)
x = coordinates[:,0]*np.cos(coordinates[:,1])
y = coordinates[:,0]*np.sin(coordinates[:,1])
X,Y = x.reshape(Npoints), y.reshape(Npoints)

data = np.where(coordinates[:,0]**2<=radius**2,data_0,0)

#fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))

#fig.suptitle(title)

#c1 = axs[0].pcolormesh(X,Y,data.reshape(Npoints), shading="auto")
#axs[0].plot(grid_center[0]*np.cos(grid_center[1]), grid_center[0]*np.sin(grid_center[1]), "rx")
#axs[0].set_title("Electric Field Data")
#cbar = fig.colorbar(c1, ax=axs[0])
#plt.show()

#exit()




# ANALYSIS

print("Starting analysis...")


# circular: p0=[A, mu_x, mu_y, sigma]
if True:

	title += " - Circular Gaussian Fit"
	popt, pcov = curve_fit(gauss_2d_wrapper(radius), (x,y), data, p0=[A, mu_x, mu_y, sigma])
	z_fit = gauss_2d_wrapper(radius)((x,y), popt[0],popt[1],popt[2],popt[3])

# elliptic: p0=[[A, mu_x, mu_y, sigma_x, sigma_y, alpha]
else:

	title += " - Elliptic Gaussian Fit"
	sigma_x, sigma_y, alpha = sigma, sigma, 0
	popt, pcov = curve_fit(gauss_elliptical_wrapper(radius), (x,y), data, p0=[A, mu_x, mu_y, sigma_x, sigma_y, alpha], bounds=([0, -np.inf, -np.inf, 0, 0, 0 ], [+np.inf, +np.inf, +np.inf, +np.inf, +np.inf, 2*np.pi ]))
	z_fit = gauss_elliptical_wrapper(radius)((x,y), popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])


z_res = z_fit-data


print("Parameters:\n",popt)
print("Cov:\n",pcov)



# POWER

data  = data.reshape(Npoints)
z_fit = z_fit.reshape(Npoints)
z_res = z_res.reshape(Npoints)


jacobian   = ( 1/np.cos(thetaphi_coordinates[:,0]) ).reshape((Npoints[0],Npoints[1]))
data_jacob = jacobian*data**2
fit_jacob  = jacobian*z_fit**2
residuals_jacob = jacobian*z_res**2

from scipy.integrate import simps

#if interp_factor<=1:
P1_u = [simps(data_jacob[i,:], np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
P_original =  simps(P1_u, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))

P2_u = [simps(fit_jacob[i,:], np.linspace(grid_lim[0],grid_lim[2],Npoints[0])) for i in range(Npoints[1])]
P_gauss =  simps(P2_u, np.linspace(grid_lim[1],grid_lim[3],Npoints[1]))

Rel_gauss = 100*P_gauss/P_original
print("Gaussian Reconstruction Relative Power:", Rel_gauss)


# PLOTS

print("Making plots...")

#title += r" - A = {:.02f}, $\mu_x$ = {:.02f}, $\mu_y$ = {:.02f}, $\sigma_x$ = {:.02f}, $\sigma_y$ = {:.02f}, $\alpha$ = {:.02f}".format(popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])

if plot_2d:

	fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))
	
	fig.suptitle(title)
	
	c1 = axs[0].pcolormesh(X, Y, 20*np.log10(abs(data)), shading="auto")
	axs[0].set_title("Electric Field Data")
	cbar = fig.colorbar(c1, ax=axs[0])

	axs[0].set_xlim(-radius,radius)
	axs[0].set_ylim(-radius,radius)
	
	c3 = axs[1].pcolormesh(X, Y, 20*np.log10(abs(z_fit)), shading="auto")
	axs[1].set_title("Fit")
	cbar = fig.colorbar(c1, ax=axs[1])
	
	c4 = axs[2].pcolormesh(X, Y, 20*np.log10(abs(z_res)), shading="auto", vmin=-40)
	axs[2].set_title("Residuals")
	cbar = fig.colorbar(c4, ax=axs[2])
	
	plt.show()


else:

	fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
	fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
	fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})
	
	ax1.plot_surface(X, Y, 20*np.log10(abs(data)))
	ax1.set_title(title + " - Original Data")
	ax2.plot_surface(X, Y, 20*np.log10(abs(z_fit)))
	ax2.set_title(title + " - Gaussian Fit")
	ax3.plot_surface(X, Y, 20*np.log10(abs(z_res)))
	ax3.set_title(title + " - Residuals")
	
	plt.show()


