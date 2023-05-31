# -----------------------------------------------------------
# contains all functions used to analyze Chandra data with 
# CIAO
#
# Functions order:
# 1. Plot Data
# 2. Compare to ACCEPT Cluster Sample
# 3. spectral property equations
#
# -----------------------------------------------------------
#
#
# Let's start with all necessary imports:
# 
#system
from __future__ import division

#numpy
import numpy as np
import numpy.polynomial.polynomial as poly

#matplotlib 
import matplotlib.pylab as plt

#astropy
from astropy.io import ascii
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.table import QTable

import pyregion



# -----------------------------------------------------------
# 1. PLOT DATA
# -----------------------------------------------------------

def pyregion_color(shape, saved_attrs):
    """ selects color for regions drawn """

    attr_list, attr_dict = saved_attrs
    attr_dict["color"] = "green"
    kwargs = pyregion.mpl_helper.properties_func_default(shape, (attr_list, attr_dict))

    return kwargs
def labels(ax, hdr):
    xname, yname = hdr['CTYPE1'], hdr['CTYPE2']
    lx = ax.coords[xname]
    ly = ax.coords[yname]
    
    if 'RA' in xname:
        xlabel = 'Right Ascension'
        
    elif 'OFFSET' in xname:
        xlabel = 'Offset'
        
    if 'DEC' in yname:
        ylabel = 'Declination'
        
    elif 'VRAD' in yname:
        ylabel = 'Velocity [%s]' % (hdr['CUNIT2'])
        
    else:
        xlabel, ylabel = xname, yname
        
    lx.set_axislabel(xlabel)
    ly.set_axislabel(ylabel)
    
    return None
def show_img(fits_image, regions=None, vmin=0.000002, vmax=0.0002, ylim=[996,1134], xlim=[1018,1155]):
    """ 
    Quickly shows fits image
    
    params:
        fits_image (str): path to fits image
        regions (str): file to regions drawn
        
    returns:
        fig, ax
    """

    #set up figure
    cmap = plt.cm.magma
    cmap.set_bad('black')

    #open data
    imgdata = fits.getdata(fits_image)
    imhdr = fits.getheader(fits_image)
    wcs = WCS(imhdr)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection=wcs)

    # custom zoom for ff.img.07-2.gz
    ax.imshow(imgdata, origin='lower', cmap='magma', 
                vmin=vmin, vmax=vmax)
    
    #ax clean up
    ax.grid(False)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    labels(ax, imhdr)

    #plot regions on top
    if regions is not None:
        #open region
        r = pyregion.open(regions).as_imagecoord(imhdr)
        patch_list, artist_list = r.get_mpl_patches_texts(pyregion_color)

        # plot each region on axis
        for p in patch_list:
            ax.add_patch(p)
        for t in artist_list:
            ax.add_artist(t)

    return fig, ax

def plot_fill_err(ax,  xdata,ydata, yerr, label, color='r', alpha=0.1, ls='none', marker='o'):
    ax.fill_between(xdata, ydata - yerr, ydata + yerr, color=color, alpha=alpha)
    ax.plot(xdata, ydata, c=color, ls=ls, marker=marker, label=label)
    return ax

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

    coeffs = poly.polyfit(ln_r, ln_xray_property, deg)

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

def dm_error_prop(variables, variable_errs):
    """ error popagation for multiplication/division
        params:
            variables (list)
            variables_err (list)
            """

    error_list = []

    for idv, variable in enumerate(variables):
        variable_err = variable_errs[idv]
        error_list.append(np.power(variable_err/variable, 2))
    
    return np.sqrt(np.sum(np.array(error_list), axis=0))

def exp_error_prop(variable, variable_err, power):
    """ error propagation for exponent"""
    return np.power(variable, power) * (power) * (variable_err/variable)

def fit_accept_clusters(data):
    """
    Create polynomial fits to ACCEPT cluster data
    return fits and errors
    """
    # Get the log temperature profile 
    ln_t = np.log(data['Tx'].value)
    ln_terr = np.log(data['Txerr'].value)

    #fit temperature profile
    r_fit, temperature_fit_unitless = fit_ln_space_polynomial(data, ln_xray_property=ln_t, deg=2)
    temperature_fit = temperature_fit_unitless * u.keV

    r_fit_err, t_fit_err = fit_ln_space_polynomial(data, ln_xray_property=ln_terr, deg=2)
    t_fit_err = t_fit_err * u.keV
    #fit pressure profile
    ln_p = np.log(data['Pitpl'].value)
    ln_perr = np.log(data['Perr'].value)
    _, pressure_fit_unitless = fit_ln_space_polynomial(data, ln_xray_property=ln_p, deg=3)
    pressure_fit = pressure_fit_unitless * u.dyne * u.cm**(-2)

    _, p_fit_err = fit_ln_space_polynomial(data, ln_xray_property=ln_perr, deg=3)
    p_fit_err = p_fit_err * u.dyne * u.cm**(-2)

    #fit electron density profile
    electron_density_fit = (pressure_fit / temperature_fit).to(u.cm**(-3))
    vars, var_errs = [pressure_fit.value, temperature_fit.value], [p_fit_err.value, t_fit_err.value]
    ed_fit_err = (electron_density_fit * dm_error_prop(vars, var_errs))
    
    #fit entropy profile
    entropy_fit = temperature_fit * np.power(electron_density_fit, -2/3)
    ed_en_sig = exp_error_prop(electron_density_fit, ed_fit_err, -2/3) 
    vars, var_errs = [temperature_fit.value, np.power(electron_density_fit, -2/3).value], [t_fit_err.value, ed_en_sig.value]
    en_fit_err = (entropy_fit * dm_error_prop(vars, var_errs) )
    
    #fit cooling time profile
    coolingFunction_fit = coolingFunction(temperature_fit)
    tcool_fit = ((3 / 2) * (1.89 * pressure_fit) / (electron_density_fit**2 / 1.07) / coolingFunction_fit) 
    ed_sig = exp_error_prop(electron_density_fit, ed_fit_err, 2) 
    cool_func_err = coolingFunction(t_fit_err)
    vars = [pressure_fit.value, np.power(electron_density_fit, 2).value, coolingFunction_fit.value]
    var_errs = [p_fit_err.value, ed_sig.value, cool_func_err.value]
    tcool_fit_err = (tcool_fit * dm_error_prop(vars, var_errs) )
    
    return r_fit, temperature_fit, pressure_fit, electron_density_fit, entropy_fit, tcool_fit, t_fit_err, p_fit_err, ed_fit_err, en_fit_err, tcool_fit_err
    
    
# -----------------------------------------------------------
# 3. spectral property equations
# -----------------------------------------------------------


def tcool(electron_density, tkev):
    #see mcdonald (2018) eqn 2
    proton_density = 0.92 * electron_density
    hydrogen_density = 0.83 * electron_density

    cooling_time = ((3/2) * ((electron_density + proton_density ) * tkev)
                    /(electron_density * hydrogen_density * cf.coolingFunction(tkev))).cgs
                    
    return cooling_time


