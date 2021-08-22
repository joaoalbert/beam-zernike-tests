import zernike_fit as zk
import zernike_class as zclass
import numpy as np


# Read from .grd

filepath = ""
radius = 0.03

grd_file = zclass.Grd_File(filepath)
grd_file.extract_data()
beam_data_file = grd_file.beams[0]

beam_data_file.generate_grid(verbose=True)
masked_data = beam_data_file.circular_mask(radius, "co")

uv_integrate(beam_data_file.co, thetas_matrix, Npoints, grid_lim)


# Reconstruct from .fits

#polar_coords = beam_data_file.polar_grid
#indices = pegar do fits
#C = fits
#Z = zk.zernike_values(beam_data_file.valid_polar_grid, self.indices, verbose)
#beam_data_zernike = np.dot(Z,C)
