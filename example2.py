#testing script for loading example zstack images and unmixing
import numpy as np
import util
import popmat
import pdb

# TODO: when inputting only one FP
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 
#filters must be inputed as a list of lists . ie [[cam1_filter1, cam2_filter],[cam1_filter1, cam2_filter]]
filts = [['TR-DFLY-F600-050'],['TR-DFLY-F521-038','TR-DFLY-F698-077']]
#cameras used for each exposure
cams = ['Andor_iXon','BSI_Prime_Express']
#FPs should be input as a lists of lists. [[exposure 1], [exposure2]] 
FPs = [['mScarlet'],['mNeonGreen','mIFP']] #these have to be exact matches as FPbase name
#excitation lines should be a list of lists. [[laser pair1], laser pair2]]
exc_lines = [['561'],['488','637']]# input paired laser lines #needs to be in quotes so can index into the csv
#input as list. %power for each laser line
laser_powers = [ .6, .7, .8] #Percentages 
#input as list, exposure time for each capture 
exposure_times = [100, 100]
beamsplit = [[1],[0,1]] #order of beamsplitter 
c_2d = util.specmix_matrix(datafolder,FPs, exc_lines, laser_powers, exposure_times,cameras = cams, filters= filts,beamsplitter = beamsplit)
