# ------------------------------------------------------------------------------------
# Fitting the Reduced GMOS Cube with PPxF (Penalized Pixel-Fitting; Cappellari (2006)) 
# ------------------------------------------------------------------------------------
# This script performs fitting of the reduced GMOS cube using the PPxF algorithm. 
# The fitting is done individually for each spaxel (spectral pixel) in the cube to 
# extract gas velocity, gas dispersion, and flux information for different emission lines.

# The fitting process is parallelized using threading, allowing for efficient computation 
# across multiple pixels. The resulting parameters and fluxes are stored in an HDF5 file 
# for further analysis!

# To run this script, in terminal type the following:
# python gmos_fit_data.py <z_galaxy> <run>
# z_galaxy = 0.336 and run is the nuber of your run-- can be anything
# ------------------------------------------------------------------------------------


from __future__ import print_function
import os

from astropy.io import fits

import numpy as np

import itertools
import h5py
import threading

import sys

import gmos_ppxf_singleG_SDSS as sg

#Change these to match your system
path_to_beads = '/Users/osaseomoruyi/Dropbox (Harvard University)/BeadsMultiwavelength/'
gmos_path = ''.join((path_to_beads, 'Analysis/gmosBeads/reduced/'))

#where to save the results
save_path = ''.join((path_to_beads, '/Analysis/gmosBeads/ppxf/fit_results/result_table/'))

z_galaxy = float(sys.argv[1])
run=str(sys.argv[2])

# Print everything
np.set_printoptions(threshold=np.inf, linewidth=np.nan)

reflines = np.genfromtxt(''.join((gmos_path, 'gmos_reflines_tex.dat')),names=True,dtype=None)


# opening raw IFU data for HE0045:
raw_file = gmos_path + 'SDSSJ1531+3414.acube.fits'
hdu_raw = fits.open(raw_file)

data_raw = hdu_raw[0].data # the shape is (3682, 31, 20) (spectra, y-coordinates, x-coordinates)
header_raw = hdu_raw[0].header

data_error = hdu_raw[1].data
header_error = hdu_raw[1].header


#eline_model from MUSE single Gaussian fits:
eline_model = fits.open(gmos_path + 'SDSSJ1531+3414.eline_model.fits')
eline_specs = eline_model[0].data

# IFU_shape=(31, 20)

# Initializing arrays
zfit_s = np.zeros(shape=(31, 20))
zfit_s[:,:] = np.nan
v_gas_s = np.zeros(shape=(31, 20))
v_gas_s[:,:] = np.nan
ev_gas_s = np.zeros(shape=(31, 20))
ev_gas_s[:,:] = np.nan
sig_gas_s = np.zeros(shape=(31, 20))
sig_gas_s[:,:] = np.nan
esig_gas_s = np.zeros(shape=(31, 20))
esig_gas_s[:,:] = np.nan
chi2_s = np.zeros(shape=(31, 20))
chi2_s[:,:] = np.nan
gas_flux_s = np.zeros(shape=(31, 20, 7))
gas_flux_s[:,:,:] = np.nan
gas_flux_err_s = np.zeros(shape=(31, 20, 7))
gas_flux_err_s[:,:,:] = np.nan
wave_s = np.zeros(shape=(31,20, 7))
wave_s[:,:,:] = np.nan
bestfit_s = np.zeros(shape=(31, 20, 7))
bestfit_s[:,:,:] = np.nan
c1fit_s = np.zeros(shape=(31, 20, 7))
c1fit_s[:,:,:] = np.nan
c2fit_s = np.zeros(shape=(31, 20, 7))
c2fit_s[:,:,:] = np.nan
raw_s = np.zeros(shape=(31, 20, 7))
raw_s[:,:,:] = np.nan


def gauss_fit(pos):
    save_image_dir = '/Users/osaseomoruyi/Dropbox (Harvard University)/BeadsMultiwavelength/Analysis/gmosBeads/ppxf/fit_results/plot_fit/run{}/'.format(run)
    os.makedirs(save_image_dir, exist_ok=True)

    pos_x, pos_y = pos[0], pos[1]
    gal_spec = data_raw[:,pos_y,pos_x]
    # cont_spec = continuum[:,pos_y,pos_x]
    err = data_error[:,pos_y,pos_x]
    eline_spec = eline_specs[:,pos_y,pos_x]

    line_spec = gal_spec #- cont_spec

    wave_spec = ((header_raw['CRVAL3'] + np.arange(0,header_raw['CDELT3']*(header_raw['NAXIS3']),header_raw['CDELT3'])))[:]

    pixel_file = "x{}_y{}".format(pos_x, pos_y)

    print('Running single gaussian fit: ')
    pp = sg.ppxf_indiv_1('SDSS1531_'+str(pos_x)+'_'+str(pos_y), z_init=z_galaxy, lambda_spec=wave_spec, galaxy_spec=line_spec,
                            err=err, header=header_raw, eline_mod=eline_spec, plotfile=pixel_file, reflines=reflines, image_save=save_image_dir)

    v_gas_s[pos_y,pos_x] = pp[0]
    ev_gas_s[pos_y,pos_x] = pp[1]

    sig_gas_s[pos_y,pos_x] = pp[2]
    esig_gas_s[pos_y,pos_x] = pp[3]

    np.put(wave_s[pos_y,pos_x], np.arange(0,len(wave_s[0,0,:])), pp[4])
    np.put(bestfit_s[pos_y,pos_x], np.arange(0,len(bestfit_s[0,0,:])), pp[5])

    chi2_s[pos_y,pos_x] = pp[6]

    np.put(gas_flux_s[pos_y,pos_x,:], np.arange(0,7), pp[7])
    np.put(gas_flux_err_s[pos_y,pos_x,:], np.arange(0,7), pp[8])

    return gas_flux_s, gas_flux_err_s, v_gas_s, ev_gas_s, sig_gas_s, esig_gas_s


xmin, xmax = 0, 20
ymin, ymax = 0, 31
INPUT = list(itertools.product(range(xmin, xmax),range(ymin, ymax)))

print('shape of IFU fitting section: ', np.shape(INPUT))

# Threading:
for arg in INPUT:
	x = threading.Thread(target=gauss_fit, args=(arg,))
	x.start()
	x.join()

    
# Writing the result to an hdf5 file:

ofile = ''.join((save_path, '{}_ppxffit.hdf5'.format(run)))
f = h5py.File(ofile, 'w')

f.attrs['name'] = 'SDSS1531'
f.attrs['flux_units'] = '10^{-16}$ erg s$^{-1}$ cm$^{-2}$'
f.attrs['velocity_units'] = 'km/s'



f.create_dataset('v_gas_s', data=v_gas_s) # gas velocity (single G fit)
f.create_dataset('ev_gas_s', data=ev_gas_s) # gas velocity (single G fit)
f.create_dataset('sig_gas_s', data=sig_gas_s) # gas sigma (single G fit)
f.create_dataset('esig_gas_s', data=esig_gas_s) # gas sigma (single G fit)
f.create_dataset('chi2_s', data=chi2_s) # chi2 of singleG fit

#################
# fluxes:
#['Hbeta' 'Halpha' '[SII]6716' '[SII]6731' '[OIII]5007_d' '[OI]6300_d' '[NII]6583_d']

f.create_dataset('Hb_flux_s', data=gas_flux_s[:,:,0]) # Hbeta flux (single G fit)
f.create_dataset('Ha_flux_s', data=gas_flux_s[:,:,1]) # Halpha flux (single G fit)
f.create_dataset('SII6717_flux_s', data=gas_flux_s[:,:,2]) # SII6717 flux (single G fit)
f.create_dataset('SII6730_flux_s', data=gas_flux_s[:,:,3]) # SII6730 flux (single G fit)
f.create_dataset('OIII5007_flux_s', data=gas_flux_s[:,:,4]*(3./4.)) #  flux (single G fit)
f.create_dataset('OIII4960_flux_s', data=gas_flux_s[:,:,4]*(1./4.)) #  flux (single G fit)
f.create_dataset('OI6363_flux_s', data=gas_flux_s[:,:,5]*(1./4.)) #  flux (single G fit)
f.create_dataset('OI6300_flux_s', data=gas_flux_s[:,:,5]*(3./4.)) #  flux (single G fit)
f.create_dataset('NII6583_flux_s', data=gas_flux_s[:,:,6]*(3./4.)) #  flux (single G fit)
f.create_dataset('NII6548_flux_s', data=gas_flux_s[:,:,6]*(1./4.)) #  flux (single G fit)


#################
# flux errors:

f.create_dataset('Hb_flux_err_s', data=gas_flux_err_s[:,:,0]) # Hbeta flux (single G fit)
f.create_dataset('Ha_flux_err_s', data=gas_flux_err_s[:,:,1]) # Halpha flux (single G fit)
f.create_dataset('SII6717_flux_err_s', data=gas_flux_err_s[:,:,2]) # SII6717 flux (single G fit)
f.create_dataset('SII6730_flux_err_s', data=gas_flux_err_s[:,:,3]) # SII6730 flux (single G fit)
f.create_dataset('OIII5007_flux_err_s', data=gas_flux_err_s[:,:,4]*(3./4.)) #  flux (single G fit)
f.create_dataset('OIII4960_flux_err_s', data=gas_flux_err_s[:,:,4]*(1./4.)) #  flux (single G fit)
f.create_dataset('OI6363_flux_err_s', data=gas_flux_err_s[:,:,5]*(1./4.)) #  flux (single G fit)
f.create_dataset('OI6300_flux_err_s', data=gas_flux_err_s[:,:,5]*(3./4.)) #  flux (single G fit)
f.create_dataset('NII6583_flux_err_s', data=gas_flux_err_s[:,:,6]*(3./4.)) #  flux (single G fit)
f.create_dataset('NII6548_flux_err_s', data=gas_flux_err_s[:,:,6]*(1./4.)) #  flux (single G fit)


# spectra 
f.create_dataset('wave_s', data=wave_s) #  wavelength
f['wave_s'].attrs['units'] = 'Angstroms'
f.create_dataset('bestfit_s', data=bestfit_s)






















