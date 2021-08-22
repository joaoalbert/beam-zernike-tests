from astropy.io import fits as pyfits
import zernike_fit as zk
import matplotlib.pyplot as plt
import numpy as np

filename = "coeffs_bingo_5_freqs.fits"

Npoints = [45,45]

with pyfits.open(filename) as hdul:

	#print(hdul.info())

	Coeffs = hdul[0].data # pyfits.getdata(filename) would work as well
	freqs = hdul[1].data["frequencies"]
	
	for f in range(len(freqs)):
	
		print("Extracting coefficients data from frequency {} GHz...".format(freqs[f]))
		C = Coeffs[f]

		beam=20*np.log10(abs(zk.beam_reconstruction(C,Npoints,verbose=True)))
		plt.pcolormesh(beam.reshape(Npoints))
		plt.show()





# Code below is deprecated.





#with pyfits.open(filename) as hdul:

	#hdul.info()
	
	#hdul[0].header
	#freqs = hdul[0].data
	#freqs = pyfits.getdata(filename) 
	# pyfits.getdata so pega as frequencias
	
	# hdul[i] -> contem a tabela de coeficientes para a freq i
	
	
	# Exemplo: dados da primeira frequencia
	
	#print(hdul[1].header)
	#ctable1    = hdul[1].data # tabela da primeira frequencia
	
	# metodo 1
	#clist1     = ctable1.field(0) # primeira coluna
	#betalist1  = ctable1.field(1) # segunda coluna
	#alphalist1 = ctable1.field(2) # terceira coluna
	
	# metodo 2
	#clist1     = ctable1.field("coefficients") # primeira coluna
	#betalist1  = ctable1.field("beta") # segunda coluna
	#alphalist1 = ctable1.field("alpha") # terceira coluna
	
	# metodo 3
	#clist1     = ctable1["coefficients"] # primeira coluna
	#betalist1  = ctable1["beta"] # segunda coluna
	#alphalist1 = ctable1["alpha"] # terceira coluna
	
	#print("")
	
	
	# Extraindo todos os dados
	
	#C_all = []
	
	#for i in range(len(freqs)):
	
		#print("Extraindo dados da frequência {} GHz...".format(freqs[i]))
		#C_all.append(hdul[i+1].data)
		
	#print(C_all)
		
	

