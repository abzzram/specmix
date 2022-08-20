# Spectral Modeling and Unmixing of Crosstalk (Specmix)
### This repo is for spectral unmixing of fluorescence intensities captured by the Dragonfly (DFLY) confocal microscope. Read through example.py for 4 channel unmixing. example2.py for 3 channel, and example3, for two channels. We employ a spectral model that calculates the theoretical senstivity constants for requested fluorophores based on DFLY properites
#### Defaults:
### Dual camera aquisition with paird laser excitation lines and 565 nm beam splitter (4 fluorophores: Blue|Green|Red|Far_Red)
####     Camera 1: Andor iXon (recieves <565 nm light from sample). Filters: TR-DFLY-F450-050 | TR-DFLY-F600-050
####    Camera 2: BSI_Prime_Express (recieves >565 nm light from sample). Filters: TR-DFLY-F521-038 | TR-DFLY-F698-077
####     Alternative emission filter and camera names can be optionally defined. Beam splitter usage can also be set to 'false'.
####        Note: when unmixing 3 channels, beamsplitter must be set to order of use per channel, see example2
### defautl cameras and filters can be changed, see example2 and example3 for usage. 
#### When using dyes, name sure dye information is in the dye data file, look for 'em' and 'ex' after the name of the dye. Make sure there is information on quantum yield as well. If there is no data for your specific dye, find it online, and add it to the data sheet. For an example of this, see 20220728_AR_multicam_test_unmix.py

### Starting new environment for specmix

#### Clone `specmix`
```
git clone https://github.com/abhineetram/specmix.git
cd specmix
```

#### Create environment
```
conda env create -f environment.yml
```

#### Activate environment
```
conda activate specmix
```

#### Try first example
```
python example.py
```
#### Notes on where AR left off 8/19/22: 
##### Branch7_5 stills needs to be edited to address Talon's last comments. Example.py file has an example for mScarlet. I found there is mScarlet can be repliaced with mStrawberry, mTangerine, or mTagRFP. However, mTagRFP seems to be the only one with decent QY/Brightness. 
#### "20220728_AR_multicam_test_unmix.py" is a script for unmixing images taken on 20020728 taken on the Dragonfly with George and Ivan. The scope was configure to image using the dual camera system with multi-aquisition. Dapi and Texas Red were imaged at the same time, and Bodipy and an empty far red filter was also imaged at the same time. The EMCCD (Andor) camera picked up the wavelength under 565nm. This data should live in a folder in the examledatasets folder, but the files were too large. The path to the files on the cluster is [Sorry, I no longer had access to the folder so I couldn't find it, it's somewhere in GB folder]. Info on the slide. Invitrogen F14781 FluoCells #2 BPAE cells Tubulin (g), Phalliodin (r), DAPI, LOT 1970997.
#### This script is not complete, we got to apply the transformation to the image, but we didn't check if it work, or if we did the math correctly.
