import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import zernike_fit as zk
import gc

#filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/bingo_standard/Retangular_MINUS1110_110_linear_x_980.0/spherical_grid.grd"
filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/bingo_standard/Retangular_30_MINUS304_linear_x_980.0/spherical_grid.grd"

output_path = "/home/joaoalb/Documents/Cosmologia/beam/results/gaussian_fits/radius_test/"

# *************** CHANGE THE TITLE ACCORDING TO FILE ****************
title_fmt = "Horn (30,-304) - 980MHz "
filename_fmt = "30_MINUS304_980" # fig


plot_2d = True # if False, will make surface plots

r_min = 0.005
r_max = 0.04
nr = 51




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

data = np.where(thetaphi_coordinates[:,0]**2<=radius**2,data_0,0)



P_originals = []
P_residuals = []
Rel_gauss = []
radii = np.linspace(r_min,r_max,nr)

for k in range(len(radii)):

	radius = radii[k]

	data = np.where(coordinates[:,0]**2<=radius**2,data_0,0)

	# ANALYSIS

	print("\nStarting analysis...")

	# circular: p0=[A, mu_x, mu_y, sigma]
	
	popt, pcov = curve_fit(gauss_2d_wrapper(radius), (x,y), data, p0=[A, mu_x, mu_y, sigma])
	z_fit = gauss_2d_wrapper(radius)((x,y), popt[0],popt[1],popt[2],popt[3])
	z_res = z_fit-data

	print("Parameters:\n",popt)
	print("Cov:\n",pcov)



	# POWER

	data  = data.reshape(Npoints)
	z_fit = z_fit.reshape(Npoints)
	z_res = z_res.reshape(Npoints)

	jacobian   = 1/(np.cos(thetaphi_coordinates[:,0])).reshape(Npoints)
	data_jacob = jacobian*data**2
	fit_jacob  = jacobian*z_fit**2
	residuals_jacob = jacobian*z_res**2

	from scipy.integrate import simps

	P1_u = [simps(data_jacob[i,:], np.linspace(grid_lims[0],grid_lims[2],Npoints[0])) for i in range(Npoints[1])]
	P_original =  simps(P1_u, np.linspace(grid_lims[1],grid_lims[3],Npoints[1]))

	P2_u = [simps(residuals_jacob[i,:], np.linspace(grid_lims[0],grid_lims[2],Npoints[0])) for i in range(Npoints[1])]
	P_res =  simps(P2_u, np.linspace(grid_lims[1],grid_lims[3],Npoints[1]))

	P_originals.append(P_original)
	P_residuals.append(P_res)
	Rel_gauss.append(100*(1-P_res/P_original))
	print("Gaussian Reconstruction Relative Power:", Rel_gauss[-1])


	# PLOTS

	print("Making plots...")

	fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(16,4))
	title    = title_fmt + " - Circular Gaussian Fit R={:.4f}".format(radius)
	filename = filename_fmt + "_r0_" + "{:04d}".format(int(radius*10**4)))
	fig.suptitle(title)
	
	c1 = axs[0].pcolormesh(X, Y, 20*np.log10(abs(data)), shading="auto")
	axs[0].set_title("Electric Field Data")
	cbar = fig.colorbar(c1, ax=axs[0])
	
	c3 = axs[1].pcolormesh(X, Y, 20*np.log10(abs(z_fit)), shading="auto")
	axs[1].set_title("Fit")
	cbar = fig.colorbar(c3, ax=axs[0])
	
	c4 = axs[2].pcolormesh(X, Y, 20*np.log10(abs(z_res)), shading="auto", vmin=-40)
	axs[2].set_title("Residuals")
	cbar = fig.colorbar(c4, ax=axs[2])

	axs[0].set_xlim(-radius,radius)
	axs[0].set_ylim(-radius,radius)
	
	plt.savefig(output_path + filename + ".png")
	
	
	arq = open(filename + ".txt", "w")
	arq.write("Radius:\n{}\n".format(radius))
	arq.write("Params:\n")
	arq.write("{} \t {} \t {} \t {} \n".format(popt[0], popt[1], popt[2], popt[3]))
	arq.write("Cov:\n")
	for line in pcov: arq.write("{} \t {} \t {} \t {} \n".format(line[0], line[1], line[2], line[3]))
	arq.write("Original Power:\n{}\n".format(P_original))
	arq.write("Residual Power:\n{}\n".format(P_res))
	arq.close()
	
	
	plt.cla()
	plt.close()
	
	
	del data
	del z_fit
	del z_res
	del jacobian
	del data_jacob
	del fit_jacob
	del residuals_jacob
	
	gc.collect()


arq = open(filename_fmt + "radius_test_{}pts".format(nr) + ".txt", "w")
arq.write("Radius range:\n({},{},{})\n".format(rmin,rmax,nr))
arq.write("Original Powers:\n")
for p in P_originals: arq.write("{}\t".format(p))
arq.write("\nResidual Powers:\n")
for p in P_residuals: arq.write("{}\t".format(p))#{}\n".format(P_residuals))
arq.close()


fig, axs = plt.subplots(ncols=1)
axs.plot(radii,Rel_gauss)
axs.set_xlabel("Radius (rad)")
axs.set_ylabel("Reconstructed Relative Power")
axs.set_title(title_fmt + " - Gaussian Fit")
plt.savefig(output_path + filename_fmt + "radius_test_{}pts.png".format(nr))


