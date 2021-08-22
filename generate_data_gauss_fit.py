import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def sinr_r2_2d(coords, A, mu_x, mu_y, sigma, radius=-1):
	'''
	A*(sin(r)/r)**2    NOT WORKING
	
	mu_x, mu_y -> not used.
	'''
	x,y = coords
	r2 = x**2 + y**2
	r02 = radius**2
	
	#output = np.zeros_like(r2) + A
	num = np.sin(r2*sigma)
	den = r2*sigma
	
	with np.errstate(divide='ignore', invalid='ignore'):
		output = np.true_divide(num,den)
		output[output == np.inf] = A
		output = np.nan_to_num(output, nan=A)
	
	# return A*(np.sin(r2*sigma)/(r2*sigma))**2 if r2!=0 else A
	
	if radius==-1: return output#np.piecewise(r2, [r2==0, r2!=0], [A, A*(np.sin(r2*sigma)/(r2*sigma))**2]) 
	else: return np.piecewise(r2, [r2<=r02, r2>r02], [output, 0])


def gauss_2d_wrapper(radius=-1):

	def gauss_2d(coords, A, mu_x, mu_y, sigma):
	
		x,y = coords
		x_0, y_0 = x-mu_x, y-mu_y
		
		r2 = x_0**2 + y_0**2
		r02 = radius**2
		
		if radius==-1: return A*np.exp( -r2 /sigma**2 /2 )
		else: return A*np.exp(-r2 /sigma**2 /2) * (r2<=r02) #np.piecewise(r2, [r2<=r02, r2>r02], [A*np.exp(-r2 /sigma**2 /2), 0])
		
	return gauss_2d
	


# Data generation

Npoints = [151,151]
x = np.linspace(-5.,5.,Npoints[0])
y = np.linspace(-5.,5.,Npoints[1])

A, mu_x, mu_y, sigma = 10, 0,0, 3
A_error, sigma_error = 1, 0.2

coords = []
for i in range(len(x)):
	for j in range(len(y)):
		coords.append([x[i], y[j]])
coords = np.array(coords)

x = coords[:,0]
y = coords[:,1]


radius = 4
model = gauss_2d_wrapper(radius) #sinr_r2_2d or gauss_2d

z_correct = model( (x,y), A, mu_x, mu_y, sigma)
z_error = z_correct + A_error * np.random.normal(0, sigma_error, z_correct.shape)
z_error = z_error * ( (x**2+y**2)<=radius**2 )
	
	
# Fit

A_est, mu_x_est, mu_y_est, sigma_est = 10,0,0,1/3 #A, mu_x, mu_y, sigma 

popt, pcov = curve_fit(model, (x,y), z_error, p0=[A, mu_x, mu_y, sigma])

print("Parametros (A, mu_x, mu_y, sigma):\n",popt)
print("Cov:\n",pcov)

z_fit = model((x,y), popt[0],popt[1],popt[2],popt[3])

z_res = z_fit-z_correct#error


# Plots

cut = True
if cut:
	# x fixo
	
	x0 = Npoints[1]//2-1
	yi = x0*Npoints[0]
	yf = (x0+1)*Npoints[0]
	
	fig = plt.figure(0)
	ax1 = fig.add_subplot(211)
	ax1.scatter(coords[yi:yf,1], z_error[yi:yf]  , c="black")
	ax1.scatter(coords[yi:yf,1], z_correct[yi:yf], c="b")
	ax1.scatter(coords[yi:yf,1], z_fit[yi:yf], c="g")
	
	ax2 = fig.add_subplot(212)
	ax2.scatter(coords[yi:yf,1], z_res[yi:yf], c="r")
	plt.show()
		

X = x.reshape(Npoints)
Y = y.reshape(Npoints)
z_error   = z_error.reshape(Npoints)
z_correct = z_correct.reshape(Npoints)
z_fit     = z_fit.reshape(Npoints)
z_res     = z_res.reshape(Npoints)

fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, figsize=(16,16), projection="3d")

#fig.suptitle(r"A = {:.02f}, $\mu_x$ = {:.02f}, $\mu_y$ = {:.02f}, $\sigma$ = {:.02f}".format(popt[0], popt[1], popt[2], popt[3]))

plot_func = pcolormesh
plot_func = plot_surface

c1 = axs[0][0].plot_func(X, Y, z_correct, shading="auto")
axs[0][0].set_title("z(x,y)")
#cbar = fig.colorbar(c1, ax=axs[0][0])

c2 = axs[0][1].plot_func(X, Y, z_error, shading="auto")
axs[0][1].set_title("z(x,y) + " + r"$\epsilon_i$")
#cbar = fig.colorbar(c2, ax=axs[0][1])

c3 = axs[1][0].plot_func(X, Y, z_fit, shading="auto")
axs[1][0].set_title("Fit")
#cbar = fig.colorbar(c3, ax=axs[1][0])

c4 = axs[1][1].plot_func(X, Y, z_res, shading="auto")
axs[1][1].set_title("Residuals ( Fit-z(x,y) )")
#cbar = fig.colorbar(c4, ax=axs[1][1])

plt.show()


#plt.subplot(

		
		

