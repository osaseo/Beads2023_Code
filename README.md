[Code, calculations, and notes for Omoruyi et al. 2023, ApJ]

DOI:

This repository houses all of the codes and calculations (contained in Jupyter Python notebooks as well as simple Python scripts) associated with our recent paper publishing Chandra, GMOS and ALMA observations of the cool core brightest cluster galaxy in SDSS 1531.

Specifically, it includes:

Reduction Codes
retrieve_[chandra/gmos/alma]_data.[py/sh] | Hopefully self-explanatory. If you run into trouble with these, you can access all of the data at each facility's public data archive!

chandra_unsharp_masks.py | Script to create unsharp masks in CIAO

chandra_clusterpyxt_params | Parameter file to do spectral fitting with clusterpyxt

gmos_reduce_data.py | Script to reduce the raw GMOS data files  (which you download using retrieve_chandra_data.py). These scripts must be run with [Py3D], in order to ensure complete reproduction of the data presented in this paper. 

gmos_pyparadise.py | Script to decouple and model the stellar and gas components of the galaxy using PYPARADISE, the Python version of the code PARADISE. Unfortunately, we do not have the exact parameter file used to decouple the data, so results may vary. 

gmos_ppxf.py | Script to decouple and model the stellar and gas components of the galaxy using PYPARADISE, the Python version of the code PARADISE. 

alma_reduce_data.py | Script to reduce the raw ASDMs (which you download using retrieve_alma_data.py) to calibrated measurement sets. These scripts must be run with [CASA] version 6.2.1.7, in order to ensure complete reproduction of the data cubes presented in this paper (although later versions of CASA will very likely also work, absent this guarantee). This script also calls alma_fluxcal_data.py. Note that running the reduction script will take roughly 12 hours on a reasonably high-end workstation with 64 GB of RAM and a 12 core processor. Machines with fewer cores and/or ram may take (much) longer. If you would like access to our reduced cubes, shoot me an e-mail (see below)!

alma_imaging/moment_maps.ipynb | These notebooks run through imaging calibrated measurement sets (generated with alma_reduce_data.py), and creates (e.g.) moment maps, a FITS cubes of the CO(3-2) data, etc. Hopefully self-explanatory.


Analysis Codes

chandra_imaging.ipynb: Notebook that creates unsharp masks
chandra_spectral_fitting.ipynb: A Jupyter notebook showing how we spectral fit and create spectral profiles

gmos_spectral_fitting.ipynb: A Jupyter notebook showing how we fit gaussians to the extracted GMOS IFU spectra

gmos_kinematic_maps.ipynb: A Jupyter notebook that makes kinematic maps and emission line diagrams from the IFU data

alma_spectral_fitting.ipynb: A Jupyter notebook showing how we fit gaussians to the extracted CO(3-2) spectra from the datacube. Note that this requires PySpecKit by Adam Ginsburg, which you can download here (requires Python 2.7, not 3.x). 

alma_pv_diagrams.ipynb: A Jupyter notebook showing how we create PV diagrams. 

alma_movie.py: A code that will make pretty movies of a MUSE or ALMA cube (see the movie above). It's been adapted from my code here.

Others in this repo are hopefully pretty self-explanatory.

As always, please feel free to email me with questions! osase.omoruyi@gmail.com