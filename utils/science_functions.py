# -----------------------------------------------------------
# contains all functions used to analyze the multi-wavelength
# data available for SDSS 1531
#
# Functions order:
# 0. ABOUT
# 1. ALMA
# 2. GMOS
# 3. HST
# 4: Other
# -----------------------------------------------------------
#
#
# Let's start with all necessary imports:
from __future__ import division
import sys
import os
import glob
from types import SimpleNamespace

#numpy
import numpy as np
import numpy.polynomial.polynomial as poly

#Astropy
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
from astropy.cosmology import LambdaCDM
import astropy.units as u
import astropy.constants as const
from astropy.convolution import convolve, Box1DKernel
from astropy.table import QTable

#Matplotlib
import matplotlib.pyplot as plt

#creating moment maps
import bettermoments as bm
from spectral_cube import SpectralCube

#spectral fitting
import pyspeckit

#seaborn
import seaborn as sns 

#pv stuff
from spectral_cube import SpectralCube

#scipy
from scipy.integrate import quad



# -----------------------------------------------------------
# 0. ABOUT
# -----------------------------------------------------------

sdss1531_dict = dict({'name': 'SDSS J1531+3414',
                    'ra': 232.7936938,
                    'dec': 34.2404172,
                    'radius': 2.5, 
                    'z': 0.335,
                    'cz': 100430.47})

# -----------------------------------------------------------
# 1. ALMA FUNCTIONS
# -----------------------------------------------------------

# NOTEBOOK 0

def channels_conversion(new_velocity_resolution, channel_number_40kms):
    """
    this function converts the velocity channel number in the 40 km/s 
    velocity cube to the corresponding channel in a cube with 
    different velocity resolution

    vres: velocity resoluton for this cube in km/s (eg. enter 20 for 20km/s)
    pchannel: channel number in 40km/s cube
    """
    reference_velocity = 40
    return channel_number_40kms * (reference_velocity / new_velocity_resolution) 

def bm_moment_maps(vcube_file, moment_type, velocity_offset_kms=0, kernel_width=3, 
                    smooth_poly_order=4, n_rms_channels=5, sigma_clip=3, thresh_mask_smooth=0, 
                    fchannel=18, lchannel=33, chatty=True):
    """Create Moment Maps using BetterMoments: https://bettermoments.readthedocs.io/en/latest/

    Parameters:
        moment_type: enter 'all' to make all moments at once or enter the specific 
                    type (eg. 'zeroth')
                
        velocity_offset_kms (int): shifts velocity zero point
        kernel_width (int): width of top hat kernel for spectral smoothing
        smooth_poly_order (int): ploynomial order for Savitzky-Golay filter
        n_rms_channels (int): number of channels used for calculating rms
        sigma_clip (int): number of sigma to be masked
        thresh_mask_smooth (int): number of pixels to spatially smooth  by
        fchannel (int): first velocity channel included
        lchannel (int): last velocity channel included
        chatty (True or False): prints what is occuring if True

    Returns:
        None. The moment maps are saved to where vcube_file is
    """
    #force into correct formats
    fchannel, lchannel = int(fchannel), int(lchannel)
    
    #load data
    data, velax = bm.load_cube(vcube_file)
    corrected_velocity_axis = velax + velocity_offset_kms 
    
    if velocity_offset_kms != 0: 
        print(f"Applied velocity offset of {velocity_offset_kms} km/s")
    
    #spectral smoothing
    if chatty is True:
        print("..spectrally smoothing the data with top hat kernel {} channels",
                " wide and applying a Savitzky-Golay filter using polynomial of",
                " order {}...".format(kernel_width, smooth_poly_order))
    
    smoothed_data = bm.smooth_data(data=data, smooth=kernel_width, polyorder=smooth_poly_order)
    
    #rms 
    if chatty is True:
        print("... estimating rms of the line free channels located in the first {} "
                "and last {} channels of the data cube".format(n_rms_channels, n_rms_channels))
    
    rms = bm.estimate_RMS(data=data, N=n_rms_channels) 
    rms_smoothed = bm.estimate_RMS(data=smoothed_data, N=n_rms_channels)
    
    if chatty is True:
        print('RMS = {:.3f} mJy/beam (original)'.format(rms * 1e3))
        print('RMS = {:.3f} mJy/beam (smoothed)'.format(rms_smoothed * 1e3))
    
    #masking
    if chatty is True:
        print("... no user masking...")
    
    user_mask = bm.get_user_mask(data=data, user_mask_path=None)
    
    #sigma clip and spatial smooth
    if chatty is True:
        print("...applying a sigma clip to all pixels with SN < {} sigma but we will first ",
                "spatially smooth the data by {} pixels...".format(sigma_clip, thresh_mask_smooth))
    
    threshold_mask = bm.get_threshold_mask(data=data,
                                           clip=(-np.inf, sigma_clip),
                                           smooth_threshold_mask=thresh_mask_smooth)

    if chatty is True:
        print("... selecting channels only between {} and {}...".format(fchannel, lchannel))
    
    #select channels to sum over
    channel_mask = bm.get_channel_mask(data=data,
                                       firstchannel=fchannel,
                                       lastchannel=lchannel)

    #combine all masks
    if chatty is True:
        print("...combining all masks...")
    mask = bm.get_combined_mask(user_mask=user_mask,
                                threshold_mask=threshold_mask,
                                channel_mask=channel_mask,
                                combine='and')

    masked_data = smoothed_data * mask
    
    
    #create moments: need moment_type: zeroth, first, second
    if chatty is True:
        print("... creating {} moment....".format(moment_type))
    #creatng all at once
    if moment_type == 'all':
        moment_types = ['zeroth', 'first', 'second', 'eighth', 'ninth']

        m0 = bm.collapse_zeroth(velax=velax, data=masked_data, rms=rms)
        m1 = bm.collapse_first(velax=velax, data=masked_data, rms=rms)
        m2 = bm.collapse_second(velax=velax, data=masked_data, rms=rms)
        m8 = bm.collapse_eighth(velax=velax, data=masked_data, rms=rms)
        m9 = bm.collapse_ninth(velax=velax, data=masked_data, rms=rms)

        #save
        bm.save_to_FITS(moments=m0, method=moment_types[0], path=vcube_file)
        bm.save_to_FITS(moments=m1, method=moment_types[1], path=vcube_file)
        bm.save_to_FITS(moments=m2, method=moment_types[2], path=vcube_file)
        bm.save_to_FITS(moments=m8, method=moment_types[3], path=vcube_file)
        bm.save_to_FITS(moments=m9, method=moment_types[4], path=vcube_file)

        #save all

    #creating a specific moment
    else:
        if moment_type == 'zeroth':
            moments = bm.collapse_zeroth(velax=velax, data=masked_data, rms=rms)
        elif moment_type == 'first':
            moments = bm.collapse_first(velax=velax, data=masked_data, rms=rms)
        elif moment_type == 'second':
            moments = bm.collapse_second(velax=velax, data=masked_data, rms=rms)
        elif moment_type == 'eighth':
            moments = bm.collapse_eighth(velax=velax, data=masked_data, rms=rms)
        elif moment_type == 'ninth':
            moments = bm.collapse_ninth(velax=velax, data=masked_data, rms=rms)
        
        #save speciic moment
        bm.save_to_FITS(moments=moments, method=moment_type, path=vcube_file)
    
    print("{} moment map created.".format(moment_type))
    return None


#NOTEBOOKS 2A AND 2B

#gas mass
def ds9_str(ra, dec, width, height, angle):
    """ creates ds9 str for ellipse in fk5 format

    Parameters:
        ra: Right Ascension
        dec: Declination
        width: diameter of ellipse
        height: full height of ellipse

    Returns:
        (int): returns ds9 str for ellipse in fk5 
    """
    #ds9 ellipse takes radii rather than diameters (astropy)
    width, height = (width.to(u.arcsec).value)/2, (height.to(u.arcsec).value)/2
    angle = 360 - angle
    
    ap_str = 'fk5; ellipse({}, {}, {}, {}, {})'.format(ra, dec, width, height, angle)
    return ap_str

def jybeam_to_jy(jy_beam, npix, cube):
    """ convert jy/beam to jy
    source: https://www.eaobservatory.org/jcmt/faq/how-can-i-convert-from-mjybeam-to-mjy/

    Parameters:
        jy_beam (float): flux in jy/beam
        npix (int): number of pixels in channel
        cube: SpectralCube

    Returns:
        (quantity): flux in jy
    """

    bmaj, bmin = cube.header['BMAJ'] * u.deg, cube.header['BMIN'] * u.deg
    
    #area of beam
    beam_area = (np.pi * bmaj.to(u.arcsec) * bmin.to(u.arcsec))/(4 * np.log(2)) 
    
    #number of pixels per beam
    num_pix_beam = cube.pixels_per_beam
    
    #total flux in each channel in Jy/beam
    fjy_beam = jy_beam * npix 
    fjy_conversion = u.beam/num_pix_beam
    return (fjy_beam * fjy_conversion).to(u.Jy) 

def mask_cube(cube, clip_value=2, rms=0.204):
    """ create a masked spectral cube

    Parameters:
        cube: SpectralCube
        clip_value (float): times sigma
        rms (float): RMS of cube in mJy

    Returns:
        cube > clipped value
    """

    rms = (rms * (u.mJy/u.beam)).to(u.Jy/u.beam)  
    mask = cube > clip_value * rms
    return cube.with_mask(mask)

def extract_cube_spectrum(cube, nrms=5, sigma_clip=3):
    """ Extract spectrum from spectral cube 

    Parameters:
        cube: SpectralCube
        nrms (int): Number of channels used to calculate rms
        sigma_clip(float): sigma clip

    Returns:
        flux, flux_err, spectral_axis

    """
    #calculate cube rms
    cube_rms = np.sqrt(np.nanvar(cube[:nrms, :, :])) * cube.unit
    #select emission > sigma
    cube = mask_cube(cube, clip_value=sigma_clip)
    #extract average cube spectrum
    cubespec = cube.mean(axis=(1, 2)) 
    #extract number of pixels in each channel
    nspec = cube.sum(axis=(1, 2))/cube.mean(axis=(1, 2))
    #convert from jy/beam to jy 
    cubespec_Jy = jybeam_to_jy(cubespec, nspec, cube) 
    #exract spectral axis
    spectral_axis = cubespec_Jy.with_spectral_unit(u.km/u.s).spectral_axis
    
    #calculate flux rms
    err_jy_beam = np.full(len(cubespec), cube_rms) * cube_rms.unit
    err_jy = jybeam_to_jy(err_jy_beam, nspec, cube) 
    flux_err = np.full(len(spectral_axis), (err_jy).to(u.mJy)) 
    #convert Nans to average err: can I do this?
    flux_err = np.nan_to_num(flux_err, copy=True, nan=np.nanmean(err_jy.value), 
                            posinf=None, neginf=None)
    
    # Ensure spectra in correct units
    flux, spectral_axis = cubespec_Jy.to(u.mJy).value, spectral_axis.to(u.km/u.s).value
    #convert NaNs to 0 for fluxes
    flux = np.nan_to_num(flux, copy=True, nan=0, posinf=None, neginf=None)
    return flux, flux_err, spectral_axis
    
def plot_spectrum(cube, ecolor='gray', guesses=None,  
                    ymin=-2, sigma=2, xminval=None, xmaxval=None, 
                    figsize=(12, 8), ymax=None, plot_spec=True):
    """ Plot and fit spectrum with pyspeckit
    Source: https://pyspeckit.readthedocs.io/en/latest/example_fromscratch.html

    Parameters:
        cube: SpectralCube
        fit_color (str): color of fitted line
        ecolor (str): color of errors
        guesses (list): input list of fit guesses
        annotate(True or False): True if you want to put fit values on plot
        direct (True or False): True if return integral of spectrum; False if
                                return integral of fit
        ymin, ymax (float): plot ymin, ymax
        sigma (float):  value for sigma clip
        xminval, xmaxval (float): plot xmax and xmin
        figsize (int, int): figure size

    Return:
        fig, ax, spectrum integral 
    """
    #extract spectrum from cube
    flux, flux_err, spectral_axis = extract_cube_spectrum(cube, sigma_clip=sigma)
    
    #Load data

    sp = pyspeckit.Spectrum(xarr=spectral_axis, data=flux, error=flux_err, 
                            xarrkwargs={'unit':'km/s'}, unit='mJy')
    fig, ax = None, None
    #plot spectrum
    if plot_spec is True:
        fig, ax = plt.subplots(figsize=figsize)
        plt.draw()

        sp.plotter(axis=ax, figure=fig, clear=False, linestyle='dashed')

        #plotting specifics
        yminval, ymaxval = ymin, 1.2 * np.nanmax(flux)
        if ymax:
            ymaxval=ymax
        if xminval is None:
            xminval, xmaxval = np.nanmin(spectral_axis), np.nanmax(spectral_axis) 
        
        sp.plotter(ymin=yminval,ymax=ymaxval, linestyle='dashed', color=ecolor, 
                    errstyle='fill', xmin=xminval, xmax=xmaxval)
        sp.plotter.label(xlabel='Velocity (km/s)',ylabel='Flux Density (mJy)') 

    return fig, ax, sp



def fit_spectrum(sp, fit_color='red',  ecolor='gray', guesses=None, annotate=False, 
                    direct=False, ymin=-2, xminval=None, xmaxval=None, 
                    figsize=(12, 8), ymax=None):
    """ Plot and fit spectrum with pyspeckit
    Source: https://pyspeckit.readthedocs.io/en/latest/example_fromscratch.html

    Parameters:
        cube: SpectralCube
        fit_color (str): color of fitted line
        ecolor (str): color of errors
        guesses (list): input list of fit guesses
        annotate(True or False): True if you want to put fit values on plot
        direct (True or False): True if return integral of spectrum; False if
                                return integral of fit
        ymin, ymax (float): plot ymin, ymax
        sigma (float):  value for sigma clip
        xminval, xmaxval (float): plot xmax and xmin
        figsize (int, int): figure size

    Return:
        fig, ax, spectrum integral 
    """

    flux, spectral_axis = sp.data, sp.xarr

    if guesses is None:
        amplitude_guess = flux.max()
        center_guess = (flux * spectral_axis).sum() / flux.sum() 
        width_guess = flux.sum() / amplitude_guess / np.sqrt(2*np.pi)
        
        guesses = [amplitude_guess, center_guess, width_guess]

    sp.specfit.multifit(guesses=guesses, negamp=False, renormalize='auto', 
                        color=fit_color)

    #plot spectrum
    fig, ax = plt.subplots(figsize=figsize)
    plt.draw()

    sp.plotter(axis=ax, figure=fig, clear=False, linestyle='dashed')
    sp.specfit.plot_components(add_baseline=False,component_yoffset=-0.2)
    sp.specfit.plotresiduals(axis=sp.plotter.axis,clear=False,yoffset=0.0,
                                label=False, color='orange')

    #plotting specifics
    yminval, ymaxval = ymin, 1.2 * np.nanmax(flux)
    if ymax:
        ymaxval=ymax
    if xminval is None:
        xminval, xmaxval = np.nanmin(spectral_axis), np.nanmax(spectral_axis) 
    
    sp.plotter(ymin=yminval,ymax=ymaxval, linestyle='dashed', color=ecolor, 
                errstyle='fill', xmin=xminval, xmax=xmaxval)
    sp.specfit.plot_fit(annotate=annotate, lw=2.0)
    sp.plotter.label(xlabel='Velocity (km/s)',ylabel='Flux Density (mJy)') 


    #extract fit results
    components = sp.specfit.get_components()

    # extract integral (units of mJy * km/s)
    integral, error = sp.specfit.integral(direct=direct, return_error=True)
    
    
    return integral, error, components

def residuals(fit, observed):
    return observed - fit
    
def molecular_mass(integral, z=0.335, xco=2e20, xco_gal=2e20):
    """ Calculates the Molecular Gas Mass from spectrum

    Parameters:
        integral (float): integral of spectrum
        z (float): redshift of galaxy
        xco (float): xco factor from Bolatto (2013)
        xco_gal (float): typically value for milky way

    Return:
        mass in solMasses
    """
    #set up cosology
    cosmo = LambdaCDM(H0=71, Om0=0.27, Ode0=0.73)
    dl = cosmo.luminosity_distance(z)
    
    #break equation up into parts
    part_a = (1.05e4/0.8)
    
    xco_units = (u.s/(u.cm * u.K * u.km))
    part_b = ((xco * xco_units)/(xco_gal * xco_units))
    
    part_c = (1/(1+z))
    
    integral = integral.to((u.Jy * u.km)/u.s)
    part_d = (integral/(u.Jy * u.km * (1/u.s))).cgs
    
    part_e = (dl/u.Mpc)**2
    
    return part_a * part_b * part_c * part_d * part_e * u.solMass
    

def mass_time(integral, integral_error, sfr=None, z=0.335, xco=2e20, 
                xco_gal=2e20, chatty=True):
    """ Calculate Molecular Hydrogen Mass and Gas Depletion time

    Parameters:
        integral (float): integral of spectrum
        integral_err (float): errors for integral of spectrum
        sfr (float): star formation rate of system
        z (float): redshift of galaxy
        xco (float): xco factor from Bolatto (2013)
        xco_gal (float): typically value for milky way
        chatty (True or False): Set True to Print values

    Return:
        H2 mass, H2 mass errro, H2 mass depletion time
    """
    
    #calculate h2 mass
    m_mol = molecular_mass(integral, z, xco, xco_gal)
    m_h2 = 0.735 * m_mol 

    #calculate h2_mass error
    m_mol_err = molecular_mass(integral_error, z, xco, xco_gal)
    m_h2_err = 0.735 * m_mol_err 

    if chatty is True:
        print('H2 mass: {:.2E} +/- {:.2E}'.format(m_h2, m_h2_err))
    
    #calculate gas depletion
    if sfr is not None:
        sfr = sfr * u.solMass/u.year
        tdep = m_h2/sfr
        
        if chatty is True:
            print(sfr)
            print('t_dep : {:.2E} '.format(tdep.to(u.Gyr)))
    else:
        tdep = None
    
    return(m_h2, m_h2_err, tdep)

def multi_gaussian_fit(xdata, mus, sigmas, amps):
    """
    Perform a multi-Gaussian fit on the given data.

    Args:
        xdata (array-like): The x-values of the data.
        mus (array-like): The means of the Gaussian components.
        sigmas (array-like): The standard deviations of the Gaussian components.
        amps (array-like): The amplitudes of the Gaussian components.

    Returns:
        array-like: The multi-Gaussian fit values.

    """

    mg_fit = gauss(xdata, mus[0], sigmas[0], amps[0]) #first gaussian

    for idm, mu in enumerate(mus):
        if idm > 0:
            mg_fit += gauss(xdata, mus[idm], sigmas[idm], amps[idm]) #add other gaussians

    return mg_fit

def multi_gauss_integral(spectral_axis, mus, sigmas, amps):
    """
    Perform the integral of a multi-Gaussian function over the given spectral axis.

    Args:
        spectral_axis (array-like): The spectral axis values.
        mus (array-like): The means of the Gaussian components.
        sigmas (array-like): The standard deviations of the Gaussian components.
        amps (array-like): The amplitudes of the Gaussian components.

    Returns:
        float: The integral of the multi-Gaussian function.

    """
    #first gaussian
    mg_integral = gauss_integral(min(spectral_axis), max(spectral_axis),  mu=mus[0], sigma=sigmas[0], amp=amps[0])[0]

    for idm, mu in enumerate(mus):
        if idm > 0:
            mg_integral += gauss_integral(min(spectral_axis), max(spectral_axis),  mu=mus[idm], sigma=sigmas[idm], amp=amps[idm])[0]
    
    return mg_integral

#Notebook 4:

def dust_mass(S_nu, nu, D, T_d):
    """
    Calculate the dust mass given the flux density, frequency, distance, and dust temperature.

    Args:
        S_nu (float or Quantity): The flux density.
        nu (float or Quantity): The frequency.
        D (float or Quantity): The distance.
        T_d (float or Quantity): The dust temperature.

    Returns:
        Quantity: The dust mass in solar masses.

    """
    kappa_nu = 0.051 * u.m**2 / u.kg

    B_nu = 2.0 * const.h * nu**3 / const.c**2 / (np.exp(const.h * nu / const.k_B / T_d) - 1.0)
    M_d = S_nu * D**2 / (kappa_nu * B_nu)
    
    return M_d.to(u.solMass)

# -----------------------------------------------------------
# 2. GMOS FUNCTIONS
# -----------------------------------------------------------

def ppxf_flux_maps(f, line_name):
    """ Plot PPXF Flux Map

    Parameters:
        f (object): file with fit
        line_name (str): Emission line name e.g for Halpha give 
                            'Ha' , for NII6583, give 'NII6583'

    Return:
        flux, flux_err
    """
    flux_s = f['{}_flux_s'.format(line_name)][()]
    flux_err_s = f['{}_flux_err_s'.format(line_name)][()]
    return flux_s, flux_err_s

def ppxf_vel_maps(f, err=False):
    """ Plot PPXF Velocity and Components Maps

    Parameters:
        f (object): file with fit

    Return:
        LOS velocity, FWHM, Chi squared
    """
    if err is False:
        vLOS = f['v_gas_s']
        chi_squared = f['chi2_s']
        fwhm = f['sig_gas_s']

    if err is True:
        vLOS = f['ev_gas_s']
        chi_squared = f['chi2_s']
        fwhm = f['esig_gas_s']
    return vLOS, fwhm, chi_squared

def ppxf_masked_flux_maps(f, line_name, snr=3.5, err=False):
    """ Mask Given Map accoridng to snr and chi square noise"""

    flux_map, flux_err_map = ppxf_flux_maps(f, line_name)
    chi_squared = f['chi2_s']

    mask = (((flux_map / flux_err_map) < snr) )
                # | (chi_squared[()] < 1) 
                # | (chi_squared[()] > 5 * np.nanmean(chi_squared[()])))
    masked_flux = np.array(flux_map).copy()
    masked_flux_err = np.array(flux_err_map).copy()

    masked_flux[mask] = np.nan
    masked_flux_err[mask] = np.nan

    if err is False:
        return masked_flux, mask
    if err is True:
        return masked_flux, masked_flux_err, mask

def ppxf_masked_vel_maps(f, mask, err=False):
    """ Apply flux map mask to given velocity maps"""

    vLOS, fwhm, chi_squared = ppxf_vel_maps(f, err=err)
    masked_vLOS = np.array(vLOS).copy()
    masked_fwhm = np.array(fwhm).copy()

    masked_vLOS[mask] = np.nan
    masked_fwhm[mask] = np.nan
    return masked_vLOS, masked_fwhm

def ebv(flux_ha, flux_hb):
    """ Calculate extinction from Ha and Hb flux. Returns E(B-V)"""

    # ka, kb = 2.63, 3.71
    # Robs = flux_ha/flux_hb
    
    # top = 2.5 * np.log10(Robs/2.86)
    # bottom = kb - ka
    
    # ebv = top/bottom

    ratio_observed = flux_ha / flux_hb
    ratio_intrinsic = 2.86
    k_alpha = 2.63
    k_beta = 3.71

    ebv = (2.5 / (k_beta - k_alpha)) * np.log10(ratio_observed / ratio_intrinsic)
    return ebv
def av(flux_ha, flux_hb):
    """ Calculate extinction from Ha and Hb flux. Returns Av"""

    e_bv = ebv(flux_ha, flux_hb)

    return 4.05 * e_bv

def lum(flux, z, cosmology):
    """ Calculate Luminosity given Flux, redshift (z),
    and defined cosmology (astropy cosmology) """
    dl = cosmology.luminosity_distance(z)
    return flux * 4 * np.pi * dl**2

def sfr_ks(flux_Ha, z, cosmology, lumin=False):
    """ Calculate SFR from Kennicut(1998) Eqn 2:
    https://iopscience.iop.org/article/10.1086/305588
    returns SFR in solar masses per year
    """
    if lumin is True:
        luminosity_Ha = flux_Ha
    if lumin is False:
        luminosity_Ha = lum(flux_Ha, z, cosmology)
    return (luminosity_Ha/(1.26e41 * u.erg/u.s)).cgs * u.solMass/u.year

def log_ne(R):
    """ Log (number of electrons) given R"""
    return 0.0543 * np.tan(-3.0553*R + 2.8506) + 6.98 - 10.6905*R + 9.9186*R**2 - 3.5442*R**3


def plot_pyparadise_results(cube_or_rss, prefix, pixel, redshift, zgas=0, mask_cont=None, 
                    fit_line=None, draw_spectral_lines=True, custom_y_limits=None, fs=16, 
                    custom_x_limits=None, fig_save_path=None, paradise_app_path=None, 
                    fsize=(13,7)):
    '''
    run like this: 
    plot_results(cube, prefix='SDSSJ1531+3414', pixel=(10,14), 
                    redshift=0.335, custom_y_limits=(-0.001,0.002), 
                    mask_cont='excl.cont', fit_line='lines.fit')
    '''
    #Load Paradise App
    sys.path.insert(0,paradise_app_path) # You'll need to fix this for your path, obviously

    try: 
        import PyParadise
        import ParadiseApp
        print('ParadiseApp successfully imported')
    except ImportError as error:
        print(f'ERROR: {error}')

    try:
        rss = PyParadise.spectrum1d.loadSpectrum(cube_or_rss)
    except IOError:
        print("Input data not found. Please check the file names.")
        sys.exit()
    try:
        cont = PyParadise.spectrum1d.loadSpectrum(f'{prefix}.cont_model.fits')
        res = PyParadise.spectrum1d.loadSpectrum(f'{prefix}.cont_res.fits')
    except IOError:
        print("PyParadise output data cannot be found. Please check the prefix for the file names.")
        sys.exit()
    try:
        line = PyParadise.spectrum1d.loadSpectrum(f'{prefix}.eline_model.fits')
        line_model_present = True
    except IOError:
        print("No emission line model found. I'll plot this without the line model.")
        line_model_present = False

    try:
        if type(pixel) == str:
            if ',' in pixel:
                cube_x = int(pixel.split(',')[0])-1
                cube_y = int(pixel.split(',')[1])-1
        elif type(pixel) == tuple:
                cube_x = int(pixel[0]-1)
                cube_y = int(pixel[1]-1)
        else:
            rss_fiber = int(pixel)-1
    except:
        print("Wrong format to select the spectrum. Please check your data format and syntax for the spectrum coordinates")
        sys.exit()

    z = redshift
    select = rss._error > 1e3
    rss._error[select] = 0

    i = 0
    plt.style.use('seaborn')
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['axes.edgecolor'] = 'k'
    
    plt.close() # just in case
    fig, (ax1, ax2) = plt.subplots(2,1,  figsize=fsize)

    if rss._datatype == 'CUBE':
        ax1.plot(rss._wave, rss._data[:, cube_y,
                 cube_x],  lw=1.5, label='data')
        ax1.plot(rss._wave, rss._error[:, cube_y,
                 cube_x],  lw=1.5, label='error')
        # if rss._mask is not None:
        #     ax1.plot(rss._wave, rss._mask[:, cube_y, cube_x]
        #              * 0.1, lw=0.5, alpha=0.8, label='badpix')
        ax1.plot(
            cont._wave, cont._data[:, cube_y, cube_x], '-r', lw=1.5, label='best-fit cont')
    elif rss._datatype == 'RSS':
        ax1.plot(rss._wave, rss._data[rss_fiber,
                 :],  lw=1.5, label='data')
        ax1.plot(rss._wave, rss._error[rss_fiber,
                 :],  lw=1.5, label='error')
        # if rss._mask is not None:
        #     ax1.plot(rss._wave, rss._mask[rss_fiber, :]
        #              * 0.1, lw=0.5, alpha=0.8, label='badpix')
        ax1.plot(cont._wave, cont._data[rss_fiber, :],
                  lw=1.5, label='best-fit cont')

    if mask_cont is not None:
        cont_mask = PyParadise.parameters.CustomMasks(mask_cont)
        for i in range(len(cont_mask['rest_frame'])):
            if i == 0:
                masklabel = mask_cont + ' (Rest frame)'
            else:
                masklabel = None
            ax1.axvspan(cont_mask['rest_frame'][i][0]*(1+z), cont_mask['rest_frame']
                        [i][1]*(1+z), color='k', alpha=0.2, label=masklabel)
        for i in range(len(cont_mask['observed_frame'])):
            if i == 0:
                masklabel = mask_cont + ' (Obs. frame)'
            else:
                masklabel = None
            ax1.axvspan(cont_mask['observed_frame'][i][0], cont_mask['observed_frame']
                        [i][1], color='r', alpha=0.2, label=masklabel)

    leg = ax1.legend(loc='upper left', fontsize=fs)
    leg.draw_frame(True)

    if line_model_present:
        if rss._datatype == 'CUBE':
            ax2.plot(res._wave, res._data[:, cube_y, cube_x],
                      lw=1.5, label='Emission Line Data (stellar subtraction residuals)')
            ax2.plot(line._wave, line._data[:, cube_y, cube_x],
                      lw=1.5, label='Emission Line Model')
        elif rss._datatype == 'RSS':
            ax2.plot(res._wave, res._data[rss_fiber, :],
                      lw=1.5, label='Emission Line Data (stellar subtraction residuals)')
            ax2.plot(line._wave, line._data[rss_fiber, :],
                      lw=1.5, label='Emission Line Model')

        if fit_line is not None:
            line_mask = PyParadise.parameters.CustomMasks(fit_line)
            for i in range(len(line_mask['rest_frame'])):
                if i == 0:
                    masklabel = fit_line + ' (Rest frame)'
                else:
                    masklabel = None
                ax2.axvspan(line_mask['rest_frame'][i][0]*(1+z), line_mask['rest_frame']
                            [i][1]*(1+z), color='b', alpha=0.2, label=masklabel)
        leg2 = ax2.legend(loc='upper left', fontsize=fs)
        leg2.draw_frame(True)
    ax1.set_ylabel(r'$f_\lambda$', fontsize=fs*1.2)
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='major', direction='in',
                    width=1.5, length=6, labelsize=16)
    ax1.tick_params(axis='both', which='minor',
                    direction='in', width=1.5, length=3)
    ax1.set_xticklabels([])
    ax2.minorticks_on()
    ax2.tick_params(axis='both', which='major', direction='in',
                    width=1.5, length=6, labelsize=16)
    ax2.tick_params(axis='both', which='minor',
                    direction='in', width=1.5, length=3)

    ax2.set_xlabel('observed-frame wavelength [$\AA$]', fontsize=fs*1.2)
    ax2.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$) ', fontsize=fs*1.2)
    
    if draw_spectral_lines is True:
        spectral_line_dict = {
                     r'H$\beta$': 4862.68,
                     r'[OIII] $\lambda 5007$': 5006.84,
                      r'H$\alpha$ + [NII]': 6562.80,
                     }
        dxval, yval = [10, 10, -250], [0.0015, 0.0005, 0.0015]
        for idl, line in enumerate(spectral_line_dict):
            ax2.axvline(spectral_line_dict[line] * (1 + zgas), lw=0.8, alpha=0.8, color='slategray')
            ax2.text(spectral_line_dict[line] * (1 + zgas) + dxval[idl], yval[idl], f'{line}')
    
    if custom_y_limits is not None:
        ax1.set_ylim(custom_y_limits)
        ax2.set_ylim(custom_y_limits)

    if custom_x_limits is not None:
        ax1.set_xlim(custom_x_limits)
        ax2.set_xlim(custom_x_limits)
    # plt.tight_layout()
    if fig_save_path is not None:
        plt.savefig(fig_save_path, dpi=300)
    plt.show()

    return fig, ax1, ax2

def ax_plot_single_gaussian_fit(ax1, save_path, pixel_file, xlim=[6450, 6750], ylim=[None, None],
                                plot_emodel=False, smooth_val=3, plot_names=False):
    """ Plots PPXF's single gaussian fit for a given pixel

    Parameters:
        save_path (str): where fit files are saved
        pixel_file (str): name of pixel file
        xlim [int, int]: where to zoom in x-axis
        ylim [int, int]: where to zoom in y-axis

    Returns:
        None
    """
    
    #Load saved files
    maskedgalaxy = np.load(save_path + '{}_sg_masked_galaxy.npy'.format(pixel_file))
    lowSN = np.load(save_path + '{}_sg_lowsn.npy'.format(pixel_file))
    wave = np.load(save_path + '{}_sg_wave.npy'.format(pixel_file))
    goodpixels = np.load(save_path + '{}_sg_goodpixels.npy'.format(pixel_file))
    pp_bestfit = np.load(save_path + '{}_sg_pp_bestfit.npy'.format(pixel_file))
    lambda_spec = np.load(save_path + '{}_sg_lambda_spec.npy'.format(pixel_file))
    zfit_gas = np.load(save_path + '{}_sg_zfitgas.npy'.format(pixel_file))
    eline_mod = np.load(save_path + '{}_sg_elinemod.npy'.format(pixel_file))
    sigma_gas = np.load(save_path + '{}_sg_sigma_gas.npy'.format(pixel_file))
    esigma_gas = np.load(save_path + '{}_sg_esigma_gas.npy'.format(pixel_file))
    noise = np.load(save_path + '{}_sg_noise.npy'.format(pixel_file))
    reflines = np.load(save_path + '{}_sg_reflines.npy'.format(pixel_file))

    gas = np.load(save_path + '{}_sg_gas.npy'.format(pixel_file))
    object_id = 'SDSS1531'

    #extract labels from each fit
    if reflines is not None:
        i_s = []
        ws = []
        labels = []

        for i,w,label in zip(range(len(reflines)),reflines['wave'],reflines['label']):
            i_s.append(i)
            ws.append(w)
            labels.append(label)
    
    #plot the fit

    #plotting smoothed spectrum
    smoothing_fact = smooth_val

    ax1.plot(wave,convolve(maskedgalaxy, Box1DKernel(smoothing_fact)),color='Gray',
                linewidth=0.5)
    ax1.plot(wave[goodpixels],convolve(maskedgalaxy[goodpixels], Box1DKernel(smoothing_fact)),
                'k', linewidth=1.)

    # overplot best fit alone
    ax1.plot(wave, pp_bestfit, 'r', linewidth=1.0,alpha=0.75, label=label)

    if plot_emodel is True:
        if len(lambda_spec) != len(eline_mod):
            ax1.plot(lambda_spec[:-1]/(1.0+float(zfit_gas)), eline_mod, 'g', linewidth=0.8, 
                        alpha=0.75, label='eline_model')
        else:
            ax1.plot(lambda_spec/(1.0+float(zfit_gas)), eline_mod, 'g', linewidth=0.8, 
                        alpha=0.75, label='eline_model')

    ax1.set_ylabel('Flux')
    ax1.set_xlabel('Rest Frame Wavelength [$\AA$]')
    #ax1.legend(loc='upper right',fontsize=10)

    ax1.set_xlim(xlim[0], xlim[1])
    ax1.set_ylim(ylim[0], ylim[1])

    #plot line names
    if plot_names is True:
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        for id_i, i in enumerate(i_s):
            w = ws[id_i]
            label = labels[id_i]

            if ((w > xmin) and (w < xmax)):
                ax1.text(w,ymin+(ymax-ymin)*(0.03+0.08*(i % 2)),'$\mathrm{'+label.decode("utf-8")+'}$',fontsize=10,\
                            horizontalalignment='center',\
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
                ax1.plot([w,w],[ymin,ymax],':k',alpha=0.5)

    return ax1

def plot_single_gaussian_fit(save_path, pixel_file, xlim=[6450, 6750], ylim=[None, None], plot_res=True):
    """ Plots PPXF's single gaussian fit for a given pixel

    Parameters:
        save_path (str): where fit files are saved
        pixel_file (str): name of pixel file
        xlim [int, int]: where to zoom in x-axis
        ylim [int, int]: where to zoom in y-axis

    Returns:
        None
    """
    
    #Load saved files
    maskedgalaxy = np.load(save_path + '{}_sg_masked_galaxy.npy'.format(pixel_file))
    lowSN = np.load(save_path + '{}_sg_lowsn.npy'.format(pixel_file))
    wave = np.load(save_path + '{}_sg_wave.npy'.format(pixel_file))
    goodpixels = np.load(save_path + '{}_sg_goodpixels.npy'.format(pixel_file))
    pp_bestfit = np.load(save_path + '{}_sg_pp_bestfit.npy'.format(pixel_file))
    lambda_spec = np.load(save_path + '{}_sg_lambda_spec.npy'.format(pixel_file))
    zfit_gas = np.load(save_path + '{}_sg_zfitgas.npy'.format(pixel_file))
    eline_mod = np.load(save_path + '{}_sg_elinemod.npy'.format(pixel_file))
    sigma_gas = np.load(save_path + '{}_sg_sigma_gas.npy'.format(pixel_file))
    esigma_gas = np.load(save_path + '{}_sg_esigma_gas.npy'.format(pixel_file))
    noise = np.load(save_path + '{}_sg_noise.npy'.format(pixel_file))
    reflines = np.load(save_path + '{}_sg_reflines.npy'.format(pixel_file))

    gas = np.load(save_path + '{}_sg_gas.npy'.format(pixel_file))
    object_id = 'SDSS1531'

    #extract labels from each fit
    if reflines is not None:
        i_s = []
        ws = []
        labels = []

        for i,w,label in zip(range(len(reflines)),reflines['wave'],reflines['label']):
            i_s.append(i)
            ws.append(w)
            labels.append(label)
    
    #plot the fit

    fig = plt.figure(figsize=(12,7))
    ax1 = fig.add_subplot(211)

    #plotting smoothed spectrum
    smoothing_fact = 3

    ax1.step(wave,convolve(maskedgalaxy, Box1DKernel(smoothing_fact)),color='Gray',
                where='mid',linewidth=0.5)
    ax1.step(wave[goodpixels],convolve(maskedgalaxy[goodpixels], Box1DKernel(smoothing_fact)),
                'k',where='mid',linewidth=1.)

    label = "Best fit template from emission lines at z={0:.3f}".format(zfit_gas)

    # overplot best fit alone
    ax1.plot(wave, pp_bestfit, 'r', linewidth=1.0,alpha=0.75, label=label)

    if len(lambda_spec) != len(eline_mod):
        ax1.plot(lambda_spec[:-1]/(1.0+float(zfit_gas)), eline_mod, 'g', linewidth=0.8, 
                    alpha=0.75, label='eline_model')
    else:
        ax1.plot(lambda_spec/(1.0+float(zfit_gas)), eline_mod, 'g', linewidth=0.8, 
                    alpha=0.75, label='eline_model')

    ax1.set_ylabel('Flux')
    ax1.set_xlabel('Rest Frame Wavelength [$\AA$]')
    ax1.legend(loc='upper right',fontsize=10)
    ax1.set_title(' '.join((object_id, pixel_file)))

    ax1.set_xlim(xlim[0], xlim[1])
    ax1.set_ylim(ylim[0], ylim[1])

    xmin, xmax = ax1.get_xlim()

    if plot_res is True:

        ax2 = fig.add_subplot(413, sharex=ax1, sharey=ax1)

        #plotting emission lines if included in the fit

        ax2.plot(wave, gas, 'b', linewidth=2, label = '$\sigma_{gas}$'+'={0:.0f}$\pm${1:.0f} km/s'.format(sigma_gas, esigma_gas))

        if len(lambda_spec) != len(eline_mod):
            ax2.plot(lambda_spec[:-1]/(1.0+float(zfit_gas)), eline_mod, 'g', linewidth=0.8, 
                        alpha=0.75, label='eline_model')
        else:
            ax2.plot(lambda_spec/(1.0+float(zfit_gas)), eline_mod, 'g', linewidth=0.8, alpha=0.75, 
                        label='eline_model')

        ymin, ymax = ax1.get_ylim()
        if (ymin < -0.5): ymin = -0.5

        ax1.set_ylim(ymin,ymax)
        ax2.set_ylim(ymin,ymax)

        ax1.set_xlim(xmin,xmax)
        ax2.set_xlim(xmin,xmax)


        ax2.set_ylabel('Best Fits')
        ax2.set_yticks(ax2.get_yticks()[::2])

        ax2.legend(loc='upper left',fontsize=10)

        # Plotting the residuals
        ax3 = fig.add_subplot(817, sharex=ax1)
        ax3.plot(wave[goodpixels], (convolve(maskedgalaxy, Box1DKernel(smoothing_fact))-pp_bestfit)[goodpixels], 
                    'k',label='Fit Residuals')
        ax3.set_yticks([-0.5,0,0.5])
        ax3.set_ylabel('Residuals')

        ax4 = fig.add_subplot(818, sharex=ax1)

        ax4.plot(wave, noise, 'k',label='Flux Error')

        ax4.set_ylabel('Noise')
        ax4.set_xlabel('Rest Frame Wavelength [$\AA$]')
        ax4.set_yticks(np.arange(0,0.5,0.1))

        ax3.set_xlim(xmin,xmax)
        ax4.set_xlim(xmin,xmax)

        for id_i, i in enumerate(i_s):
            w = ws[id_i]
            label = labels[id_i]

            if ((w > xmin) and (w < xmax)):
                ax1.text(w,ymin+(ymax-ymin)*(0.03+0.08*(i % 2)),'$\mathrm{'+label.decode("utf-8")+'}$',fontsize=10,\
                            horizontalalignment='center',\
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
                ax1.plot([w,w],[ymin,ymax],':k',alpha=0.5)

        fig.subplots_adjust(hspace=0.03)
    return None

def Classify_nii_BPT(table):
    """ Classify NII BPT spaxels"""
    snr = 10

    SNR_Ha, Ha_mask = ppxf_masked_flux_maps(table, 'Ha', snr=snr)
    SNR_Hb, Hb_mask = ppxf_masked_flux_maps(table, 'Hb', snr=snr)
    SNR_OIII, O3_mask = ppxf_masked_flux_maps(table, 'OIII5007', snr=snr)
    SNR_NII, N2_mask = ppxf_masked_flux_maps(table, 'NII6583', snr=snr)


    ## BPT DIAGRAM: NII ##
    #Kewley et al. 2001: starburst vs AGN classification. Solid lines in BPT
    #log10(flux_oiii_5006/flux_hbeta)=0.61/(log10(flux_nii_6583/flux_halpha)-0.47)+1.19
    #Kauffmann et al. 2003: starburst vs composites. Dashed line in BPT
    #log10(flux_oiii_5006/flux_hbeta)=0.61/(log10(flux_nii_6583/flux_halpha)-0.05)+1.3
    #Schawinsky et al. 2007: Seyferts vs LINERS
    #log10(flux_oiii_5006/flux_hbeta)=1.05*log10(flux_nii_6583/flux_halpha)+0.45

    i_bptnii = np.log10(SNR_NII/SNR_Ha)
    j_bptnii = np.log10(SNR_OIII/SNR_Hb)
    Kew01_nii = 0.61/(i_bptnii-0.47)+1.19
    Scha07 = 1.05*i_bptnii+0.45
    Ka03 = 0.61/(i_bptnii-0.05)+1.3
    
    agn_nii = ((j_bptnii>=Kew01_nii) & (j_bptnii>Scha07))
    liner_nii = ((j_bptnii>=Kew01_nii) & (j_bptnii<Scha07) | (i_bptnii>=0.47))
    composite_nii = ((j_bptnii>=Ka03) & (j_bptnii<Kew01_nii))
    sf_nii = ((j_bptnii<Ka03) & (i_bptnii<=-0.25))

    excitation = np.zeros_like(SNR_Ha)
    excitation[sf_nii] = 1
    excitation[composite_nii] = 2
    excitation[liner_nii] = 3
    excitation[agn_nii] = 4
    excitation[excitation==0.0] = np.nan

    
    return (i_bptnii, j_bptnii, excitation)

# generate Gaussian function
def gauss(x, mu, sigma, amp):
    return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def gauss_integral(lowerbound, upperbound, mu, sigma, amp):
    integral = quad(gauss, lowerbound, upperbound, args=(mu, sigma, amp))
    return integral

# -----------------------------------------------------------
# 3. HST FUNCTIONS
# -----------------------------------------------------------

def ysc_load(path_to_beads):
    """ Load Young Stellar Superclusters Table"""
    table_path = path_to_beads + 'Analysis/tables/'
    cluster_radec = np.load(table_path + 'ysc_coords.npy')
    return cluster_radec

def load_HST_data(path_to_beads, id_file=0, fprint=False):
    """ Loads HST data
    Paramteres:
        id_file (int): 0- f606
                        1 - f814
                        2 - f390
                        3 - f160
    """

    #HST 
    hst_files = glob.glob(''.join((path_to_beads, '/Analysis/hst_data/hst_*')))

    ff = hst_files[id_file]
    if fprint is True:
        print(ff)

    fhdu = fits.open(ff)[0]

    hst_hdr = fhdu.header
    hst_wcs = WCS(hst_hdr)
    hst_hdu = fhdu
    return hst_hdr, hst_wcs, hst_hdu


# -----------------------------------------------------------
# 4. OTHER
# -----------------------------------------------------------

# -----------------------------------------------------------
# 4A. OTHER ALMA FUNCTIONS
# -----------------------------------------------------------

#OTHER:

def rho(mass, a, r):
    #equation 1
    top = mass * a
    bottom = 2 * np.pi * np.power((r + a), 3)
    return top/bottom

def scale_length(r_half):
    #equation 2
    factor = 2 + np.sqrt(2)
    return r_half/factor

def potential(a, r, mass):
    #equation 3
    top = -1 * const.G * mass
    bottom = r + a
    return top/bottom

def v_freefalll(mass, r, a, r0=10*u.kpc, v0=0):
    #equation 4
    term1 = np.power(v0,2)/2
    term2 = potential(a, r0, mass)
    term3 = potential(a, r, mass)
    return np.sqrt(2 * (term1 + term2 - term3))

# def scale_length(r_half):
#     #equation 2
#     #factor = 2 + np.sqrt(2)
#     return r_half/1.8153#r_half/factor

def distance(angular_size, hdr, cosmology, zh=0.335, hdr_field=None):
    """calculate the  angular distance, given angular size.

        Parameters:
            angular_size: with any astropy units
            hdr: regular fits header
            cosmology (object): defined astropy cosmology 
        Returns:
            (quantitiy) distance in kpc
    """
    angular_size = angular_size.to(u.deg)
    
    da = cosmology.angular_diameter_distance(zh)
    dl = cosmology.luminosity_distance(zh)
    
    if hdr_field:
        pix_size = hdr[hdr_field]* u.deg
    else:
        pix_size = hdr['CDELT2']* u.deg

    rpix = ((pix_size * da) /(206265 * u.arcsec)).to(u.kpc).value * u.kpc
    return (rpix * (angular_size/pix_size)).to(u.kpc)

def freefall_velocity(height, hdr, cosmology, halo_mass=1e13*u.solMass, 
                        stellar_mass=2*np.power(10, 11.58)*u.solMass,
                        num_intervals=25, v0=610*(u.km/u.s), radius=25*u.kpc,
                        radius_0=0*u.kpc):
    """ 
    calculates velocity and distance travelend given a height
    
    """
    halo_mass = 1e13 * u.solMass #estimate within 30kpc using figure 3
    stellar_mass = np.power(10, 11.58) * u.solMass * 2 #stellar mass of one elliptical * 2

    r_effective = radius/2
    a = scale_length(r_effective) #guess half mass radius is half of the radius

    radial_distances = np.linspace(0, radius.value, num_intervals) * u.kpc
    radial_vels = v_freefalll(stellar_mass, radial_distances, a, r0=height, 
                                v0=v0).to(u.km/u.s)- v0

    #make all negative radial velocities 0
    idr = np.where((radial_vels > 0) & (~np.isnan(radial_vels)))
    radial_vels = radial_vels[idr]
    radial_distances = radial_distances[idr] + radius_0


    #convert velocity and distance to pixels so can plot
    wcs = WCS(hdr)
    pixel_vels = [int(wcs.wcs_world2pix(0, 0, radial_vel * 1000, 0)[2]) for radial_vel in radial_vels]
    pixel_dist = radial_distances/distance(hdr['CDELT2'] * u.deg, hdr, cosmology)

    return pixel_vels, pixel_dist

def vel_cube(scube_file, vel_convention='radio'):
    """This function converts a spectralcube to a velocity cube"""
    #create velocity cube 
    cube = SpectralCube.read(scube_file)
    vcube = cube.with_spectral_unit(u.km / u.s, velocity_convention=vel_convention)
    return vcube


def distance(angular_size, hdr, cosmology, zh=0.335, hdr_field=None):
    """calculate the  angular distance, given angular size.

        Parameters:
            angular_size: with any astropy units
            hdr: regular fits header
            cosmology (object): defined astropy cosmology 
        Returns:
            (quantitiy) distance in kpc
    """
    angular_size = angular_size.to(u.deg)
    
    da = cosmology.angular_diameter_distance(zh)
    dl = cosmology.luminosity_distance(zh)
    
    if hdr_field:
        pix_size = hdr[hdr_field]* u.deg
    else:
        pix_size = hdr['CDELT2']* u.deg

    rpix = ((pix_size * da) /(206265 * u.arcsec)).to(u.kpc).value * u.kpc
    return (rpix * (angular_size/pix_size)).to(u.kpc)

# -----------------------------------------------------------
# 4B. OTHER GMOS FUNCTIONS
# -----------------------------------------------------------

def fov_info(fov_file):
    hdu = fits.open(fov_file) # Open the .fits file, return a list of Header/Data Units (HDUs). 
    
    hdr = hdu[0].header #The header of the data fits file and it tells you all technical info.
    fovdata = hdu[0].data # literally the data
    dim = hdu[0].data.shape[1:] #telling us the shape of the data think (x,y)
    
    w = WCS(hdr)

    return fovdata, hdr, dim, w

#read emission line table
def gas_info(gas_file):
    # Read in the emission line table, from which we'll make gas maps
    eline_hdu = fits.open(gas_file) #opening the data
    
    eline_tab = eline_hdu[1].data
    eline_columns = eline_hdu[1].header
    
    print(eline_columns)
    
    eline_x_cor = eline_tab.field('x_cor')
    eline_y_cor = eline_tab.field('y_cor')
    
    return eline_tab, eline_x_cor, eline_y_cor

#creating stellar kinematics maps
def kinematics_maps(stars_file, dim, cz, subtract=None):
    
    """
    stars_file: filename of stellar hdu
    dim: dimensions from fov cube
    cz= speed of light * z
    """
    #read in table
    stellar_hdu = fits.open(stars_file)
    stellar_tab = stellar_hdu[1].data #trying to pick out the data
    stellar_columns = stellar_hdu[1].header # selecting the header
    
    print("stellar table shape: {}".format(stellar_tab.shape))
    
    # Make arrays for the spatial coordinates, velocity and fwhm
    stellar_x_cor, stellar_y_cor = stellar_tab.field('x_cor'), stellar_tab.field('y_cor') 
    stellar_vel, stellar_vel_err = stellar_tab.field('vel_fit'), stellar_tab.field('vel_fit_err')
    stellar_fwhm, stellar_fwhm_err = stellar_tab.field('disp_fit'), stellar_tab.field('disp_fit_err')
    
    xdim_stars = np.max(stellar_x_cor) - np.min(stellar_y_cor)
    ydim_stars = np.max(stellar_y_cor) - np.min(stellar_y_cor)
    dim_stars = (xdim_stars, ydim_stars)
    print("xdim * ydim = {}".format(xdim_stars * ydim_stars))
    
    # mask the data, selecting S/N > 500 spaxels. 
    stellar_select = ((stellar_vel / stellar_vel_err) > 5.0) 
    
    #match data to fov dimensions and make empty maps of NaNs
    stellar_vel_map = np.full((dim[0], dim[1]), np.nan)
    stellar_disp_map = np.full((dim[0], dim[1]), np.nan)

    print("Map shape: {}".format(np.shape(stellar_vel_map)))
    
    #make stellar maps
    stellar_vel_map[stellar_y_cor[stellar_select], stellar_x_cor[stellar_select]] = stellar_vel[stellar_select]
    systemic = cz.to(u.km/u.s).value 
    print("Median vel: {}".format(np.nanmedian(stellar_vel_map)))
    
    if subtract:
        stellar_vsys_map = stellar_vel_map
        
    else:
        stellar_vsys_map = stellar_vel_map - systemic #systemic velocity subtracted stellar velocities
        
    #make fwhm maps
    # We're in *such* a high S/N regime that I'm just going to use the same stellar_select 
    # boolean mask as above, even though I technically created it based on S/N in the *velocity*
    # fit, not with FWHM fit. The Velocity and FWHM maps will therefore have an identical footprint. 
    # Really, everything is so high S/N that this is basically not relevant. 
    print("created stellar and fwhm maps")
    stellar_fwhm_map = np.full((dim[0], dim[1]), np.nan)
    stellar_fwhm_map[stellar_y_cor[stellar_select], stellar_x_cor[stellar_select]] = stellar_fwhm[stellar_select]
    
    maps = stellar_vsys_map, stellar_fwhm_map
    
    return maps

#Creating GMOS emission line maps
def emission_maps(line_name, err, eline_tab, eline_y_cor, eline_x_cor, dim, cz):
    """
    To create diff maps:
    name - filament name e.g 'Hbeta_flux' --> just give 'Hbeta' 
    err - desired S/N ratio
    dim - fov cube dimensions
    """

    #flux map
    flux = eline_tab.field(line_name + '_flux')
    
    flux_map = np.full((dim[0],dim[1]) ,np.nan)
    if err==None:
         
        flux_map[eline_y_cor,eline_x_cor] = flux
    else:
        #this is a scientific way to cut out noise.
        flux_err = eline_tab.field(line_name + '_flux_err')
        gas_select = ((flux / flux_err) > err)
        
        flux_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = flux[gas_select]
    
    
    #vel map
    vel = eline_tab.field(line_name + '_vel')
    if cz is None:
        median_vel = np.nanmedian(vel)
    else:
        median_vel = cz.value
    
    vel_map = np.full((dim[0],dim[1]) ,np.nan)
    
    if err==None:
        vel_map[eline_y_cor,eline_x_cor] = vel - median_vel
    
    else:
        vel_err = eline_tab.field(line_name + '_vel_err')
        gas_select = ((vel / vel_err) > err)

        vel_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = vel[gas_select] - median_vel

    #fwhm map
    fwhm = eline_tab.field(line_name + '_fwhm')
    fwhm_map = np.full((dim[0],dim[1]), np.nan)
    
    if err==None:
        fwhm_map[eline_y_cor,eline_x_cor] = fwhm

    else:
        fwhm_err = eline_tab.field(line_name + '_fwhm_err')
        gas_select = ((fwhm / fwhm_err) > err)

        fwhm_map[eline_y_cor[gas_select],eline_x_cor[gas_select]] = fwhm[gas_select]
    #print(np.nanmax(fwhm_map.flatten()))

    return [flux_map, vel_map, fwhm_map]


#----
# CHANDRA FUNCTIONS
#----



# -----------------------------------------------------------
# 2. COMPARE TO ACCEPT Cluster Sample
# The Archive of Chandra CLuster Entropy Propfile (ACCEPT) is a 
# database storing key physical properties of observed galaxy 
# clusters derived from the X-ray data. 
# 
# This is just for underplotting profiles against our SDSS 1531 
# profiles (i.e. for error checking and comparison). The below 
# functions fit an n-degree polynomial to the ACCEPT main table 
# and computes relevant quantities. Ask Grant if you have questions. 
# This is explained in more detail [in this 
# repository](https://github.com/granttremblay/RainMaker). 
# https://web.pa.msu.edu/astro/MC2/accept/
# -----------------------------------------------------------

def parse_accept_table(accept_table_file):
    """ read ACCEPT Table and asssign units"""

    # Read in the raw (unitless) accept table as an astropy table
    accept_table_raw = ascii.read(accept_table_file)

    # assign units to each column
    Name = accept_table_raw['Name']
    Rin = accept_table_raw['Rin'] * u.Mpc
    Rout = accept_table_raw['Rout'] * u.Mpc
    nelec = accept_table_raw['nelec'] * u.cm**(-3)
    neerr = accept_table_raw['neerr'] * u.cm**(-3)
    Kitpl = accept_table_raw['Kitpl'] * u.keV * u.cm**2
    Kflat = accept_table_raw['Kflat'] * u.keV * u.cm**2
    Kerr = accept_table_raw['Kerr'] * u.keV * u.cm**2
    Pitpl = accept_table_raw['Pitpl'] * u.dyne * u.cm**(-2)
    Perr = accept_table_raw['Perr'] * u.dyne * u.cm**(-2)
    Mgrav = accept_table_raw['Mgrav'] * u.M_sun
    Merr = accept_table_raw['Merr'] * u.M_sun
    Tx = accept_table_raw['Tx'] * u.keV
    Txerr = accept_table_raw['Txerr'] * u.keV
    Lambda = accept_table_raw['Lambda'] * u.erg * u.cm**3 / u.s
    tcool52 = accept_table_raw['tcool5/2'] * u.Gyr
    tcool52err = accept_table_raw['t52err'] * u.Gyr
    tcool32 = accept_table_raw['tcool3/2'] * u.Gyr
    tcool32err = accept_table_raw['t32err'] * u.Gyr

    names = ('Name', 'Rin', 'Rout', 'nelec', 'neerr', 'Kitpl',
             'Kflat', 'Kerr', 'Pitpl', 'Perr', 'Mgrav', 'Merr',
             'Tx', 'Txerr', 'Lambda', 'tcool52', 't52err',
             'tcool32', 't32err'
             )

    # Create the new table as a QTable with units 
    main_accept_table = QTable(
        [Name, Rin, Rout, nelec, neerr, Kitpl,
         Kflat, Kerr, Pitpl, Perr, Mgrav, Merr,
         Tx, Txerr, Lambda, tcool52, tcool52err,
         tcool32, tcool32err], names=names
    )

    return main_accept_table



def fit_ln_space_polynomial(data, ln_xray_property, deg):
    """Fit a DEG-order polynomial in extrapolated x, y space.

    numpy.polynomial.polinomial.polyfit() returns coefficients,
    from 0th order first to N-th order last (note that this is
    *opposite* from how np.polyfit behaves!).

    Params: 
        data (QTable): table from parse_accept_table
        ln_xray_property: natural log of the given property to 
                        fit (eg. Temperature)
        deg: degrees of freedom in polynomial fit to data (eg. 2)
    """

    # Obtain radii
    r = (data['Rin'] + data['Rout']) * 0.5
    ln_r = np.log(r.value)

    # Generate  wider ange of radii in log space
    log10_r_range = np.arange(300) / 100 - 3

    # un-log10 radii, give it a unit
    r_range = (10**log10_r_range) * u.Mpc

    # Also give its unitless natural log, used for fitting
    # with polyval() and fit_polynomial()'s coefficients
    ln_r_range = np.log(r_range.value)

    idf = np.isfinite(ln_r) & np.isfinite(ln_xray_property)
    coeffs = poly.polyfit(ln_r[idf], ln_xray_property[idf], deg)

    # polyval() is used to assemble cubic fit:
    # $p(x) = c_0 + c_1 x + c_2 x^2 + c3 x^3$
    # where c_n are the coeffs returned by polyfit()
    ln_fit = poly.polyval(ln_r, coeffs)
    fit = np.exp(ln_fit)

    # Now use these coefficients to extrapolate fit
    # across larger radius

    ln_fit_range = poly.polyval(ln_r_range, coeffs)
    fit_range = np.exp(ln_fit_range)

    fitpackage = (fit, r, fit_range, r_range, coeffs)

    return r_range, fit_range


def coolingFunction(kT):
    r"""
    Implement the Tozzi & Norman (2001) cooling function.

    This is an analytic fit to Sutherland & Dopita (1993), shown
    in Equation 16 of Parrish, Quataert, & Sharma (2009),
    as well as Guo & Oh (2014).

    See here: arXiv:0706.1274. The equation is:

    $\Lambda(T) = [C_1 \left( \frac{k_B T}{\mathrm{keV}} \right)^{-1.7}
                  + C_2\left( \frac{k_B T}{\mathrm{keV}} \right)^{0.5}
                  + C_3] \times 10^{-22}$
    """
    # For a metallicity of Z = 0.3 Z_solar,
    C1 = 8.6e-3
    C2 = 5.8e-2
    C3 = 6.3e-2

    alpha = -1.7
    beta = 0.5

    # kT.value is supposed to be unitless,
    # the kT in the real equation is divided by keV. I'm not cheating!
    coolingFunction = (C1 * kT.value**alpha
                       + C2 * kT.value**beta
                       + C3
                       ) * 1.0e-22 * (u.erg * u.cm**3 / u.s)

    return coolingFunction

def dm_error_prop(variables, variable_errs, function):
    """ error popagation for multiplication/division
        params:
            variables (list)
            variables_err (list)
            """

    error_list = []

    for idv, variable in enumerate(variables):
        variable_err = variable_errs[idv]
        error_list.append(np.power(variable_err/variable, 2))
    
    return function * np.sqrt(np.sum(np.array(error_list), axis=0))

def exp_error_prop(variable, variable_err, power):
    """ error propagation for exponent"""
    return np.power(variable, power) * (power) * (variable_err/variable)


def M_r(T, rho_g, r):
    """ 
    Vikhlinin (2006) eqn 7
    """
    dln_rho_g = np.diff(np.log(rho_g)) / np.diff(np.log(r)) 

    dln_T = np.diff(np.log(T)) / np.diff(np.log(r))

    return -3.68e13 * T[:-1] * r[:-1] * (dln_rho_g + dln_T) 

def t_ff(mass, radius):
    return np.sqrt((2 * radius**3)/(const.G * mass))