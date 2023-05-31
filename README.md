[Code, calculations, and notes for Omoruyi et al. 2023, ApJ]

DOI:

This repository houses all of the codes and calculations (contained in Jupyter Python notebooks as well as simple Python scripts) associated with our recent paper publishing Chandra, GMOS and ALMA observations of the cool core brightest cluster galaxy in SDSS 1531.

Specifically, it includes:

**Reduction Codes: data_reduction/**

*[chandra/gmos/alma]_retrieve_data.[py/sh]* | Hopefully self-explanatory. If you run into trouble with these, you can access all of the data at each facility's public data archive!


*chandra_unsharp_masks.py* | Script to create unsharp masks in CIAO

*chandra_clusterpyxt_params* | Parameter file to do spectral fitting with clusterpyxt

*gmos_fit_data.py* | Script using PPxF to fit the gas component of the galaxy created from PYPARADISE, the Python version of the code PARADISE. This script makes use of gmos_ppxf_singleG_SDSS.py and gmos_reflines_tex.dat

*alma_reduce_imaging/run_all_cleans.py* | Script to reduce the raw ASDMs (which you download using retrieve_alma_data.py) to calibrated measurement sets. These scripts must be run with [CASA] version 6.2.1.7, in order to ensure complete reproduction of the data cubes presented in this paper (although later versions of CASA will very likely also work, absent this guarantee). Note that running the reduction script will take roughly 12 hours on a reasonably high-end workstation with 64 GB of RAM and a 12 core processor. Machines with fewer cores and/or ram may take (much) longer. If you would like access to our reduced cubes, send me an e-mail (see below)!

*alma_best_imaging/moment_maps.ipynb* | These notebooks run through imaging calibrated measurement sets (generated with alma_reduce_data.py), and creates (e.g.) moment maps, a FITS cubes of the CO(3-2) data, etc. Hopefully self-explanatory.


**Analysis Codes: data_analysis/**

Each of the scripts/notebooks below make use of functions found in the **/utils** folder

*chandra/*

*chandra_imaging.ipynb* | Notebook that visualizes surface brightness maps and unsharp masks created fromt he reduction stage 

*chandra_spectral_profiles.ipynb* | Notebook that shows how we create spectral profiles in CIAO

*chandra_spectral_maps.ipynb* | Notebook that plots the results from ClusterPyxt

*chandra_analysis.ipynb* | Notebook that analyzes the chandra data further

*gmos/*

*gmos_eline_analysis.ipynb* | Notebook that analyzes the optical emission lines from the fit GMOS IFU spectra

*alma/*

*alma_gas_mass.ipynb* | Notebook showing how we fit gaussians to the extracted CO(3-2) spectra from the datacube. Note that this requires PySpecKit by Adam Ginsburg. 

*alma_poisition_velocity_diagrams.ipynb* | Notebook showing how we create PV diagrams. 

*alma_analysis.py* | Notebook that further analyzes the ALMA data. 

*hst/*

*hst_analysis.ipynb* | Notebook that analyzes the HST data further. 


*multi_wavelength_comparison.ipynb* | Notebook that compares the multi-wavelength data.

Please feel free to email me with any questions! osase.omoruyi@gmail.com