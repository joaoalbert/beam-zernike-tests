import zernike_class as zclass
import zernike_fit as zk
import numpy as np



def txt_extract_coeffs(arq):
	'''
	simplified version
	'''
	
	brackets = 3
	coeffs = []
	
	line = arq.readline()
	while "]"*brackets+"\n" not in line:
		line = arq.readline()
		data = line.replace("[","").replace("]","").split()
		coeffs.append([ float(d) for d in data ])
		
	return np.array([coeffs])
	
	

#filepath = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/bingo_standard/Retangular_30_MINUS304_linear_x_980.0/"
filepath = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/asymetric/Retangular_MINUS930_MINUS304_linear_x_980.0/"

# Fast
#filepath = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/npoints_test/980_151pts_04uv/"


#results_path = "/home/joaoalb/Documents/Cosmologia/beam/results/correct_grid/horn_MINUS1110_110/uv005_pts1001_r003_b22/"
#results_path = "/home/joaoalb/Documents/Cosmologia/beam/results/correct_grid/horn_30_MINUS304/GZ_analysis_r003_b22"

results_path = "/home/joaoalb/Documents/Cosmologia/beam/results/tests/"



verbose = True


zernike_run=False
if zernike_run:

	radius   = 0.1 # 0.025
	beta_max = 22
	pol = "co"

	indices    = None
	show_plot  = True
	title      = "freq980"
	fig_path   = results_path + title + ".png"
	final_txt  = results_path + title + ".txt"
	final_fits = results_path + title + ".fits"

	zernike_analysis = zclass.Zernike_Analysis(radius, beta_max, filepath=filepath, verbose=verbose, pol=pol, indices=indices, show_plot=show_plot, fig_path=fig_path, final_txt=final_txt, final_fits=final_fits)
	
	
	
gauss_fit=True
if gauss_fit:

	grd_file = zclass.Grd_File(filepath)
	grd_file.extract_data(verbose)
	beam_data = grd_file.beams[0]
	
	show_plot  = True
	title      = "gauss_freq980"
	fig_path   = results_path + title + ".png"
	final_txt  = results_path + title + ".txt"
	
	radius = 0.03
	
	beam_gaussian_fit = zclass.Beam_Gaussian_Fit(beam_data, radius, verbose, show_plot, fig_path, final_txt)
	
	
	
	
GZ_analysis=False
if GZ_analysis:

	radius   = 0.03 # 0.025
	beta_max = 22
	pol = "co"

	indices    = None
	show_plot  = True
	gauss_title   = "gauss_freq980"
	zernike_title = "zernike_freq980"
	gauss_fig_path     = results_path + gauss_title + ".png"
	zernike_fig_path   = results_path + zernike_title + ".png"
	final_txt = results_path + ".txt"
	#zernike_final_fits = results_path + zernike_title + ".fits"

	zclass.Gaussian_Zernike_Analysis(filepath, radius, beta_max, verbose=verbose, show_plot=show_plot, gauss_fig_path=gauss_fig_path, zernike_fig_path=zernike_fig_path, final_txt=final_txt, pol=pol) #final_fits=final_fits
	
	
	
	
	
recalculate_beam=False
if recalculate_beam:

	txt_file = "/home/joaoalb/Documents/Cosmologia/beam/results/correct_grid/horn_30_MINUS304/GZ_analysis_r003_b22/GZ_analysis_r003_b22.txt"
	txt = open(txt_file, "r")	
	for i in range(31): txt.readline()
	coeffs = txt_extract_coeffs(txt)[0]
	txt.close()
	
	Npoints = np.array([301, 301])
	grid_lims = np.array([-0.0548, -0.0027, 0.0452, 0.0973])
	grid_center = np.array([-0.0044, 0.0475])
	cols = np.zeros((Npoints[0]*Npoints[1],4))
	frequency = 0.98
	
	beam_data = zclass.Beam_Data(cols, grid_lims, grid_center, Npoints, frequency)
	beam_data.generate_grid(verbose=True, fast=True)
	masked = beam_data.circular_mask(0.03, "co")
	
	coords = beam_data.valid_polar_grid
	Z = zk.zernike_values(coords, coeffs[:,1:3], verbose=True)
	valid_beam = np.dot(Z,coeffs[:,0])
	beam = zk.reconstruct_masked(valid_beam, masked).reshape(Npoints)

	
	import matplotlib.pyplot as plt
	
	X = (beam_data.polar_grid[:,0]*np.cos(beam_data.polar_grid[:,1])).reshape(Npoints)
	Y = (beam_data.polar_grid[:,0]*np.sin(beam_data.polar_grid[:,1])).reshape(Npoints)
		
	plt.pcolormesh(X,Y,20*np.log10(abs(beam)),shading="auto", vmin=-90, vmax=30)
	#plt.xlim(X[0,0],X[-1,-1])
	#plt.ylim(Y[0,0],Y[-1,-1])
	plt.colorbar()
	plt.show()
	
	
	
