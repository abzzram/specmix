from calendar import c
import numpy as np
import util
import popmat
import pdb

#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 

# FPs = [['mTagBFP2','TagRFP-T'],['mNeonGreen','iRFP670']] #these have to be exact matches as FPbase name
FPs = [['mTagBFP2','mScarlet'],['mNeonGreen','iRFP670']] #these have to be exact matches as FPbase name

#alternatives to scarlet are mStrawberry, mTangerine, but QY is low
exc_lines = [['405','561'],['488','637']]# input paired laser lines #needs to be in quotes so can index into the csv

laser_powers = [1,1,1,1] #Percentages 
exposure_times = [100, 200]
c_2d = util.specmix_matrix(datafolder,FPs, exc_lines, laser_powers, exposure_times)

#noramlize matrix for printing
row_sums = c_2d .sum(axis=1)
new_matrix = c_2d  / row_sums[:, np.newaxis]

print('\n Channel-fluorophore cross talk matrix (Normalized): \n')
np.set_printoptions(suppress=2, precision=3)
print(new_matrix)

# pdb.set_trace()
abundances = np.array([1, 0.5, 0.2, 0.7])
intensities = c_2d @ abundances
#add noise to intensities
noise = np.random.normal(0, 10, intensities.shape)
noisy_intensities = intensities + noise
noisy_intensities
# inverse
abundance_est = np.linalg.pinv(c_2d) @ noisy_intensities

print('simulated abundances:')
print(abundances)
print('adding noise to simulated intensities...')
print('\n')
print('predicted abundances from unmixing:')
print(abundance_est)

