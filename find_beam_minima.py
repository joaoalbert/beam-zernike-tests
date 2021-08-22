import numpy as np
import matplotlib.pyplot as plt
import zernike_fit as zk



# asymmetric
#filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/bingo_standard/Retangular_MINUS1110_110_linear_x_980.0/spherical_grid.grd"
# symmetric
filename = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/bingo_standard/Retangular_30_MINUS304_linear_x_980.0/spherical_grid.grd"

cols, grid_lims, grid_centers, Nps, freqs = zk.grasp_spher_grd_file(filename, shift_center=True, verbose=False)	

grid_lim = grid_lims[0]
grid_center = grid_centers[0]
Npoints = Nps[0]
cols = cols[0]

data_0 = np.sqrt(cols[:,0]**2 + cols[:,1]**2)
log_data = 20*np.log10(abs(data_0))
data_0 = data_0.reshape(Npoints)
i_max, j_max = np.unravel_index(data_0.argmax(), data_0.shape)



# filtering those weird structures
# that appear because of the horizontal/vertical minima-search pattern
filtering = False
fraction = 25
filter_range_hor = list(range(i_max - Npoints[0]//fraction, i_max + Npoints[0]//fraction + 1))
filter_range_ver = list(range(j_max - Npoints[1]//fraction, j_max + Npoints[1]//fraction + 1))
if not filtering:
	filter_range_hor = []
	filter_range_ver = []



from scipy.signal import find_peaks

#horizontal selection
minimum_indices_hor   = [ find_peaks(-1*data_0[i])[0] for i in range(Npoints[0])]
minimum_indices_hor   = [ [idx for idx in idxs if idx not in filter_range_hor ] for idxs in minimum_indices_hor ]
minimum_indices_hor_u = [ minimum_indices_hor[i][j] for i in range(len(minimum_indices_hor)) for j in range(len(minimum_indices_hor[i])) ]
minimum_indices_hor_v = [ i for i in range(Npoints[0]) for j in range(len(minimum_indices_hor[i])) ]

#vertical selection
minimum_indices_ver   = [ find_peaks(-1*data_0[:,i])[0] for i in range(Npoints[1])]
minimum_indices_ver   = [ [idx for idx in idxs if idx not in filter_range_hor ] for idxs in minimum_indices_hor ]
minimum_indices_ver_v = [ minimum_indices_ver[i][j] for i in range(len(minimum_indices_ver)) for j in range(len(minimum_indices_ver[i])) ]
minimum_indices_ver_u = [ i for i in range(Npoints[1]) for j in range(len(minimum_indices_ver[i])) ]


#watershed-like algorithm

ws_min = 5.
ws_min = 10.**(ws_min/20.)
ws_below_indices   = np.arange(len(data_0.flatten()))[data_0.flatten()<=ws_min]
ws_below_indices_u = ws_below_indices%Npoints[0]
ws_below_indices_v = ws_below_indices//Npoints[0]


# pixel distances (number of pixels)
dists_1 = [ np.sqrt((i-i_max)**2 + (j-j_max)**2) for i,j in zip(minimum_indices_hor_u, minimum_indices_hor_v)  ]
dists_1.sort()
dists_2 = [ np.sqrt((i-i_max)**2 + (j-j_max)**2) for i,j in zip(minimum_indices_ver_u, minimum_indices_ver_v)  ]
dists_2.sort()
dists_3 = [ np.sqrt((i-i_max)**2 + (j-j_max)**2) for i,j in zip(ws_below_indices_u, ws_below_indices_v)  ]
dists_3.sort()



plt.figure(0)
plt.pcolormesh( 20*np.log10(abs(data_0)) )
plt.plot(minimum_indices_hor_u, minimum_indices_hor_v, "bo")
plt.plot(minimum_indices_ver_u, minimum_indices_ver_v, "ro")
plt.colorbar()

#plt.figure(1)
#plt.plot(minimum_indices_hor_u, minimum_indices_hor_v, "bo")
#plt.figure(2)
#plt.plot(minimum_indices_ver_u, minimum_indices_ver_v, "ro")

plt.figure(3)
plt.plot(dists_1,"bo")
plt.plot(dists_2,"ro")

plt.figure(4)
plt.plot(dists_3,"go")

plt.figure(5)
plt.pcolormesh( 20*np.log10(abs(data_0)) )
plt.plot(ws_below_indices_u, ws_below_indices_v, "go")
plt.colorbar()

plt.show()
exit()

#horizontal slice

row = 70 #Npoints[0]//2

plt.figure(10)
plt.plot(minimum_indices_hor[row], 20*np.log10(abs(data_0[row][minimum_indices_hor[row]])), "x")
plt.plot(20*np.log10(abs(data_0[row])))
plt.show()



