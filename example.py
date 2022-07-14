import numpy as np
import util
import popmat

#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 

paths = util.get_filepaths(datafolder)

FPs = ['mTagBFP2','mNeonGreen','mScarlet','iRFP670'] #these have to be exact matches as FPbase name
exc_lines = ['488','561']#,'514','685'] #needs to be in quotes so can index into the csv
em_filters = ['>530', '<530'] # TODO
specdata = util.get_spectra(FPs, paths, exc_lines)

# laser_powers = [10, 10, 10, 10]
# exposure_times = [100, 100]

# assemble unmixing matrix
c_2d = popmat.populate_matrix(specdata, exc_lines, em_filters, FPs)

# simulate
abundances = np.array([1, 0.5, 0.2, 0.7])
intensities = c_2d @ abundances

# inverse
abundance_est = np.linalg.pinv(c_2d) @ intensities

