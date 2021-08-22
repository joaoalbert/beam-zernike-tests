import zernike_class as zclass
import numpy as np


grd_file_path = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/alex/Retangular_MINUS150_110_linear_x_1120.0/"


# Extracting data from .grd files

grd_file = zclass.Grd_File(grd_file_path)
grd_file.extract_data(verbose=False)


# Since one .grd file may contain multiple frequencies,
# the data below are given inside lists, although we
# usually produce .grd's with only one frequency per file.

freqs        = grd_file.frequencies
grid_lims    = grd_file.grid_lims
grid_centers = grd_file.grid_centers
Nps          = grd_file.Nps

beams = grd_file.beams


# Extracting data for the first (and probably only) frequency

freq        = freqs[0]
grid_lim    = grid_lims[0]
grid_center = grid_centers[0]
Npoints     = Nps[0]

beam = beams[0]
beam_co = beam.co
beam_cx = beam.cx


# Generating the grid points

beam.generate_grid(fast=True) # This fast option must be False only for proper polar coordinates around maximum
uv_grid       = beam.uv_grid
thetaphi_grid = beam.thetaphi_grid


# Usage example

beam_matrix = beam_co.reshape(Npoints)

r   = thetaphi_grid[:,0]
phi = thetaphi_grid[:,1]

XX = (r*np.cos(phi)).reshape(Npoints)
YY = (r*np.sin(phi)).reshape(Npoints)


import matplotlib.pyplot as plt

plt.figure(1)
plt.pcolormesh(XX,YY,20*np.log10(abs(beam_matrix)), shading="auto")
plt.colorbar()
plt.show()
