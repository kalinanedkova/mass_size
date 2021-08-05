# This python code can be used to interpolate the Chebyshev polynomials from Galapagos-2
# Inputs: (i)   Chebyshev polynomial coefficients
#         (ii)  Effective wavelengths of the filters in the correct order
#         (iii) The wavelength at which the polynomial will be evaluated
# NOTE:   (1)   Please make sure entries (ii) and (iii) have the same units. 

########################################################
###################### EXAMPLES ########################
########################################################
# For the HFF Structural catalogues,  we've used 7 bands (F125W, F435W, F606W, 
# F814W, F105W, F140W, F160W) with effective wavelengths: 1.2486,  0.43273, 0.59218, 
# 0.80593, 1.0552, 1.3923,  1.5369 in microns
# Suppose you want to interpolate the effective radius Chebyshev polynomial in order to
# obtain the 4000A rest-frame size. For an example galaxy at z = 1, that means
# evaluating the Chebyshev polynomial at 0.8microns.
# Suppose RE_GALFIT_CHEB for the object in question is: 
# (5.71, -1.02, 0.28, 0.0, 0.0, 0.0, 0.0)

################# Running in a script #################
# To use this code in a script:
# import get_Cheb_poly
# get_Cheb_poly.interpolate_cheb(cheb=[5.71, -1.02, 0.28], bands=[1.2486,  0.43273,
# 0.59218, 0.80593, 1.0552, 1.3923,  1.5369], wave=0.8)
# 
# This returns a float, which is the Chebyshev polynomial value at the desired
# wavelength -- i.e., for this example: 5.83420872
#
# To interpolate entire an column of Chebyshev polynomials from a Galapagos-2 catalogue:
# import get_Cheb_poly
# get_Cheb_poly.interpolate_cheb_from_col(cheb=cat['RE_GALFIT_BAND'], bands=[1.2486,
# 0.43273, 0.59218, 0.80593, 1.0552, 1.3923,  1.5369], wave=0.4*(cat['z']+1))
# returns an array of length equal to the length of the catalogue, 'cat'
# for objects for which the polynomial is undefined, the code returns -99
# For example: array([ -99., -99., ..., 0.97521114, 2.86086448, ..., -99., -99.])

############ Running from the command line ############
# For quick calculations for single objects, you may also wish to run this from the
# command line. To do this, type:
# python get_Cheb_poly.py --cheb '5.71, -1.02, 0.28' --bands '1.2486,  0.43273, 0.59218,
# 0.80593, 1.0552, 1.3923,  1.5369' --wave '0.8'
# (For this example, the code evaluates a Chebyshev polynomial with coefficients
# 5.71, -1.02, 0.28 (note this the same as entering 5.71, -1.02, 0.28, 0.0, 0.0, 0.0, 0.0)
# This will print : At wavelength 0.8 : 5.834208721774708
# i.e., the size of this object would be ~5.8pix at 0.8microns 
# The effective wavelengths of the filters should be entered in the order which is 
# used by Galapagos-2. For the HFF structural catalogue, that is F125W, F435W, F606W, 
# F814W, F105W, F140W, and F160W

########################################################
# Note: this code can be used for Chebyshev polynomials of any order 
#       (e.g., 1st order polynomial, etc.).
#       The order is defined by the number of elements given in --cheb
# Note: you may use any number of bands, so long as they correspond to the number of 
#       bands used with Galapagos-2
# Note: Please check that the wavelength at which the polynomial is evaluated is within
#       the wavelength range of 'bands'. Otherwise, the Chebyshev polynomial is being
#       extrapolated and the results are not as reliable. The code will print a warning
#       when this happens. To ignore this warning, use: 
#       import warnings; warnings.filterwarnings("ignore")
########################################################

# numpy.polynomial.chebyshev uses Chebyshev polynomials of the first kind
from numpy.polynomial.chebyshev import Chebyshev
import numpy as np
import argparse
import warnings

def interpolate_cheb(cheb,bands,wave):
	bluest = np.min(bands) #effective wavelength of bluest filter
	reddest = np.max(bands) #effective wavelength of reddest filter
	if wave[i] >= reddest or wave[i] <= bluest:
		warnings.warn("The wavelength at which you're interpolating the Chebyshev polynomial is outside the wavelength range of the used bands")

	func = Chebyshev(cheb, domain=[bluest, reddest]) # This evaluates the Chebyshev polynomial
	return func(wave)

# To do the above for more than one object:
def interpolate_cheb_from_col(cheb_column,bands,wave):
	bluest = np.min(bands) #effective wavelength of bluest filter
	reddest = np.max(bands) #effective wavelength of reddest filter
	if len(cheb_column) != len(wave):
		warnings.warn("wave must be an array with the same number of elements as cheb_column")
	else:
		funcs = {n: Chebyshev(cheb_column[n], domain=[bluest, reddest]) for n in range(len(cheb_column))}
	
		polynom = np.zeros(len(cheb_column))
		for i in range (len(cheb_column)):
			if cheb_column[i][0] != -99: 
				polynom[i] = funcs[i](wave[i])
				if wave[i] >= reddest or wave[i] <= bluest:
					warnings.warn("The wavelength at which you're interpolating the Chebyshev polynomial for index:"+str(i)+" is outside the wavelength range of the used bands")
			else : polynom[i] = -99		
			
	return polynom


if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Interpolate Chebyshev polynomials of the first kind')
	parser.add_argument('--cheb', dest='cheb',help='input should be comma separated chebyshev polynomial coefficients')
	parser.add_argument('--wavelength', dest='wave',help='input should be the wavelength value at which the polynomial will be evaluated')
	parser.add_argument('--bands', dest='bands',help='input should be comma separated wavelengths of the used bands')

	args = parser.parse_args()

	cheb = list(map(float,args.cheb.split(',')))
	cheb = [ float(i) for i in cheb]

	bands = list(map(float,args.bands.split(',')))
	bands = [ float(i) for i in bands]

	wave = list(map(float,args.wave))
	wave = [ float(i) for i in wave]


	bluest = np.min(bands) #effective wavelength of bluest filter
	reddest = np.max(bands) #effective wavelength of reddest filter

	func = Chebyshev(cheb, domain=[bluest, reddest])
	polynom = np.zeros(len(wave))
	for i in range(len(wave)):
		polynom[i] = func(wave[i]) # This evaluates the Chebyshev polynomial
		print('At wavelength ' + str(wave[i]) + ' : ' +str(polynom[i]) +'\n')

	
