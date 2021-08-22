import numpy as np
import matplotlib.pyplot as plt


filepath = "/home/joaoalb/Documents/Cosmologia/GRASP_Student/NewProject/Job_uv0501_elliptical/spherical_grid.grd"

grd_file = open(filepath, "r")

line = grd_file.readline()
while line!="++++\n":	line = grd_file.readline()
grd_file.readline(), grd_file.readline()

# general file information
grid_center = [ float(d) for d in grd_file.readline().split() ]
grid_lim    = [ float(d) for d in grd_file.readline().split() ]
Npoints     = [ int(d)   for d in grd_file.readline().split()[:2] ]




# Now the file is divided into the "slices" of the ellipsis.
# The division lines, which are lines containing only two numbers,
# indicate the position of the first pixel and the number of 
# valid pixels in that slice.



# Here, I will fill in the blank spaces, reconstructing the
# rectangular structure. An alternative would be to simply
# save each slice information (division lines) and deal with
# the non-rectangular structure later.

masked_point = [0,0,0,0] # this will be masked later
data = []


while True:

	slice_info  = grd_file.readline().split()
	if slice_info==[]: break
	
	slice_start  = int(slice_info[0])
	slice_size   = int(slice_info[1])
	l_blank_size = slice_start-1
	r_blank_size = Npoints[0]-l_blank_size-slice_size

	# filling in the left-side blank space
	for i in range(l_blank_size): data.append(masked_point)
	
	for l in range(slice_size):
		line_data = [ float(d) for d in grd_file.readline().split() ]
		data.append(line_data)
		
	# filling in the right-side blank space
	for j in range(r_blank_size): data.append(masked_point)
		
		
data = np.ma.masked_where(data==0,data)
co = np.sqrt(data[:,0]**2 + data[:,1]**2).reshape(Npoints)

plt.pcolormesh(20*np.log10(co))
plt.show()
		
