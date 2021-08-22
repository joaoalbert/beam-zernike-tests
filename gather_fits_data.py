import astropy.io.fits as pyfits
import zernike_fit as zk
import numpy as np


filepath = "/home/joaoalb/Documents/Cosmologia/beam/results/correct_grid/horn_MINUS1110_110/uv005_pts1001_r003_b21/"
fits_fmt = filepath + "mhz{}_r003_b21_uv005_nps1001.fits"
freqs    = np.linspace(980,1260,11)[1:-1]
filelist = np.array([fits_fmt.format(int(f)) for f in freqs])



for i in range(len(filelist)):
	
	filename = filelist[i]
	
	with pyfits.open(filename) as arq:
	
		if i==0:
			coeffs = arq[0].data
			freqs = arq[1].data["frequencies"]	
		
		else:
			coeffs = np.concatenate((coeffs,arq[0].data))
			freqs  = np.concatenate((freqs,arq[1].data["frequencies"]))
			
			
C = zk.spectral_plot(coeffs, freqs, msc=2.5, verbose=True, fig_path=filepath+"spectral_representation_msc2_5.png")
