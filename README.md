# description, installation instructions, brief usage (or detailed example what dyes and FPs are available)
# This repo is for spectral unmixing of fluorescence intensities captured by the Dragonfly (DFLY) confocal microscope
# Read through example.py for 4 channel unmixing. example2.py for 3 channel, and example3, for two channels
# We employ a spectral model that calculates the theoretical senstivity constants for requested fluorophores based on DFLY properites
# Defaults:
## Dual camera aquisition with paird laser excitation lines and 565 nm beam splitter (4 fluorophores: Blue|Green|Red|Far_Red)
###     Camera 1: Andor iXon (recieves <565 nm light from sample). Filters: TR-DFLY-F450-050 | TR-DFLY-F600-050
###     Camera 2: BSI_Prime_Express (recieves >565 nm light from sample). Filters: TR-DFLY-F521-038 | TR-DFLY-F698-077
###     Alternative emission filter and camera names can be optionally defined. Beam splitter usage can also be set to 'false'.
###         Note: when unmixing 3 channels, beamsplitter must be set to order of use per channel, see example2

# Starting new environment for specmix

# Clone `specmix`
```
git clone https://github.com/abhineetram/specmix.git
cd specmix
```

# Create environment
conda env create -f environment.yml

# Activate environment
conda activate specmix

# Try first example
python example.py
