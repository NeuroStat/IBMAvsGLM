# README

Initial scripts and attempt of using resting state as null data comes from https://github.com/wanderine/ParametricMultisubjectfMRI.

# Data
For each subject, we have a 3D structural scan (with and without skull) and one 4D functional (resting state) image from the 1000 functional connectomes project.

Scanning/pre-processing parameters:
* TR(s) = 3
* nscan = 119
* Functional image: 72 x 72 x 47 with 3 x 3 x 3 mm voxels
* Structure image (without skull): 144 x 192 x 192 with 1.2 x 1.98 x 1.98 mm voxels.
* Highpass filtering = yes
* Lowpass filtering = no
* Slice timing correction = no
* Prewhitening = yes
* Motion parameters in model = yes
* Linear registration to subject main structural image = yes (applying transformation = at second level)
  * BBR cost function
  * Full search space
* Linear registration to standard space = yes (applying transformation = at second level)
  * MNI152 T1 2mm brain
  * BBR cost function
  * full search space
  * 12 DOF


# Analyses
At the moment, we run analyses with these settings:

### First level settings
* One scanning site (Cambridge)
* Block design task:
  * 10 seconds ON/OFF
* Smoothing at (mm):
  * 4
  * 6
  * 8
  * 10
* Highpass filter cut-off at 20 seconds (twice length of blocked condition)
* Square waveform for EV
* Gamma convolution
* Cluster RFT inference at level 0.05 with Z-threshold (cluster forming) at 2.3

### Second level GLM settings

> Not run at the moment.