from __future__ import print_function

# from time import clock
import numpy as np

from ppxf import ppxf
import ppxf.ppxf_util as util

# Print everything
np.set_printoptions(threshold=np.inf, linewidth=np.nan)

path_to_beads = '/Users/osaseomoruyi/Dropbox (Harvard University)/BeadsMultiwavelength/'
ppxf_path = ''.join((path_to_beads, 'Analysis/gmosBeads/ppxf/'))


reflines = np.genfromtxt(''.join((ppxf_path, 'run_fit/reflines_tex.dat')),names=True,dtype=None)


def vac_to_air(wave):
    air = wave / (1.0 + 2.735182E-4 + 131.4182 / (wave**2.) + 2.76249E8 / (wave**4.))
    return air


    #######################################################################################
###    The function below takes an individual galaxy in the lega-c survey and uses ppxf  ###
###    to fit a combination of Conroy+ SSP templates and emission lines to a lega-c ###
###    spectrum.                                                                    ###
###                                                                                 ###
#######################################################################################
##
####################################      Inputs:       ########################################
##
## object_id = the id number of the galaxy
##
## z_init = The initial guess for the spectroscopic redshift of the object. This has to be a 
##          a reasonably good guess since pPXF does a poor job of redshift estimation. Not a 
##          problem here since these are known.
##
## galaxy_spec = The galaxy flux array to be used
## 
## err = The error array to be used
##
## Optional Inputs:
##
## plotfile = The filename of an output plot of the best-fit
##
## outfile = An open file object to which to write a line including the best-fit parameters
##
## outfile_spec = An open file object to which to output the log-binned spectrum and best-fit
##                models
##
## mpoly = The order of the multiplicative polynomial (I recommend just using the defaults)
##
## apoly = The order of the additive polynomial (I recommend just using the defaults)
## 
## reflines = A catalog of reference lines to lable on the output plot
##
##
####################################      Outputs:       ########################################
##
## zfit_stars, ezfit_stars = redshift (and error) of the continuum model 
##
## sigma_stars, esigma_stars = sigma (and error) of the continuum model
##
## zfit_gas, ezfit_gas = redshift (and error) of the emission lines
##
## sigma_gas, esigma_gas = sigma (and error) of the emission lines
##
## SN_median = Median Signal-to-noise in the full spectrum
## 
## SN_rf_4000 = Signal-to-noise at rest-frame 4000 Angstroms
## 
## SN_obs_8030 = Signal-to-noise at observed-frame 8030 Angstroms
## 
## chi2 = chi^2 value from PPXF associated with the best-fit model
##

 
def ppxf_indiv_1(object_id, z_init, lambda_spec, galaxy_spec, err, header, eline_mod,
              plotfile=None, outfile=None, outfile_spec=None, mpoly=None, apoly=None,
              reflines=None, image_save=None):
    """ This functions fits one pixel only. See fitting details above."""

    # Speed of light
    c = 299792.458

    # Crop to finite
    use = np.isfinite(galaxy_spec) & (galaxy_spec > 0.0)

    # Making a "use" vector
    use_indices = np.arange(galaxy_spec.shape[0])[use]
    galaxy_spec = (galaxy_spec[use_indices.min():(use_indices.max()+1)])[2:-3]
    err = (err[use_indices.min():(use_indices.max()+1)])[2:-3]
    lambda_spec = (lambda_spec[use_indices.min():(use_indices.max()+1)])[2:-3]
    eline_mod = (eline_mod[use_indices.min():(use_indices.max()+1)])[2:-3]
    
    lamRange_gal = np.array([np.min(lambda_spec),np.max(lambda_spec)])
    FWHM_gal = 2.3548*(np.max(lambda_spec) - np.min(lambda_spec))/len(lambda_spec)

    lamRange_gal = lamRange_gal/(1+float(z_init)) # Compute approximate restframe wavelength range
    FWHM_gal = FWHM_gal/(1+float(z_init))   # Adjust resolution in Angstrom

    # log rebinning for the fits
    galaxy, logLam_gal, velscale = util.log_rebin(np.around(lamRange_gal,decimals=3), galaxy_spec, flux=True)

    noise, logLam_noise, velscale = util.log_rebin(np.around(lamRange_gal,decimals=3), err, \
                                              velscale=velscale, flux=True)

    # correcting for infinite or zero noise
    noise[np.logical_or((noise == 0.0),np.isnan(noise))] = 1.0
    galaxy[np.logical_or((galaxy < 0.0),np.isnan(galaxy))] = 0.0

    if galaxy.shape != noise.shape:
        galaxy = galaxy[:-1]
        logLam_gal = logLam_gal[:-1]
        

    # Construct a set of Gaussian emission line templates.
    # Estimate the wavelength fitted range in the rest frame.
    gas_templates, line_names, line_wave = util.emission_lines(logLam_gal, lamRange_gal, FWHM_gal)
    
    goodpixels = np.arange(galaxy.shape[0])
    wave_ = np.exp(logLam_gal)
    wave = vac_to_air(wave_)
    
    # exclude lines at the edges (where things can go very very wrong :))
    include_lines = np.where((line_wave > (wave.min()+10.0)) & (line_wave < (wave.max()-10.0)))
    if line_wave[include_lines[0]].shape[0] < line_wave.shape[0]:
        line_wave = line_wave[include_lines[0]]
        line_names = line_names[include_lines[0]]
        gas_templates = gas_templates[:,include_lines[0]]

    
    #set uo variables for ppxf fit

    reg_dim = gas_templates.shape[1:]
    dv = 0  # km/s
    templates = np.hstack([gas_templates])
    vel = 0. # Initial estimate of the galaxy velocity in km/s

    #lines and components to fit
    nNLines = gas_templates.shape[1]
    component = [0]*nNLines
    gas_component = [True]*len(component)

    start_gas = [vel, 100.]
    moments = [2] # fit (V,sig,h3,h4) for the stars and (V,sig) for the gas' broad and narrow components
    fixed = None
    start = start_gas # adopt the same starting value for both gas (BLs and NLs) and stars

    ## Additive polynomial degree
    
    if apoly is None: degree = int(np.ceil((lamRange_gal[1]-lamRange_gal[0])/(100.0*(1.+float(z_init)))))
    else: degree = int(apoly)

    if mpoly is None: mdegree = 3
    else: mdegree = int(mpoly)
    
    print('ADDITIVE POLYNOMIAL = {0} MULTIPLICATIVE POLYNOMIAL = {1}'.format(degree,mdegree))    

    bounds_gas = [[-2000,2000], [-2000,2000]]
    bounds = bounds_gas

    pp = ppxf.ppxf(templates, galaxy, noise, velscale, start, fixed=fixed,
              plot=False, moments=moments, mdegree=mdegree,
              degree=degree, vsyst=dv, reg_dim=reg_dim, gas_component=gas_component,
              component=component, goodpixels=goodpixels, bounds=bounds)

    vel_gas = pp.sol[0]
    evel_gas = pp.error[0]*np.sqrt(pp.chi2)

    sigma_gas = pp.sol[1] #first # is template, second # is the moment
    esigma_gas = pp.error[1]*np.sqrt(pp.chi2)

    zfit_gas = (z_init + 1)*(1 + pp.sol[0]/c) - 1

    ezfit_gas = (z_init + 1)*pp.error[0]*np.sqrt(pp.chi2)/c

    include_gas = True

    #for saving files to plopt the fits later
    if plotfile is not None:
        save_path = image_save
        print('Saving parameters to {0}'.format(plotfile))

        maskedgalaxy = np.copy(galaxy)
        lowSN = np.where(noise > (0.9*np.max(noise)))
        maskedgalaxy[lowSN] = np.nan

        wave = np.exp(logLam_gal)*(1.+z_init)/(1.+zfit_gas)
        # NOTE: if there's an offset, its probably due to z_init not being = to zfit_gas or something
        np.save(save_path + '{}_sg_masked_galaxy'.format(plotfile), maskedgalaxy)
        np.save(save_path + '{}_sg_lowsn'.format(plotfile), lowSN)
        np.save(save_path + '{}_sg_wave'.format(plotfile), wave)
        np.save(save_path + '{}_sg_goodpixels'.format(plotfile), goodpixels)
        np.save(save_path + '{}_sg_pp_bestfit'.format(plotfile), pp.bestfit)
        np.save(save_path + '{}_sg_lambda_spec'.format(plotfile), lambda_spec)
        np.save(save_path + '{}_sg_zfitgas'.format(plotfile), zfit_gas)
        np.save(save_path + '{}_sg_ezfitgas'.format(plotfile), ezfit_gas)
        np.save(save_path + '{}_sg_elinemod'.format(plotfile), eline_mod)
        np.save(save_path + '{}_sg_sigma_gas'.format(plotfile), sigma_gas)
        np.save(save_path + '{}_sg_esigma_gas'.format(plotfile), esigma_gas)
        np.save(save_path + '{}_sg_noise'.format(plotfile), noise)
        np.save(save_path + '{}_sg_reflines'.format(plotfile), reflines)
        # for plotting emission lines if included in fit
        if include_gas:
            gas = pp.matrix[:,-nNLines:].dot(pp.weights[-nNLines:])
            np.save(save_path + '{}_sg_gas'.format(plotfile), gas)
            
    return vel_gas, evel_gas, sigma_gas, esigma_gas, wave, pp.bestfit, pp.chi2, pp.gas_flux, pp.gas_flux_error