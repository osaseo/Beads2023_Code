# -----------------------------------------------------------
# contains all functions used to plot the multi-wavelength
# data available for SDSS 1531
#
# Functions order:
# 1. Plot Style
# 2. FIGURE COMPONENTS
# 3. FITS Handling
# 4: PLOT FITS MAPS
# 5: MOVIES
# -----------------------------------------------------------
#
#
# Let's start with all necessary imports:
from __future__ import division
import os
import glob

#numpy
import numpy as np

#matplotlib 
import matplotlib.pylab as plt
from matplotlib import ticker
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes

#astropy
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import Angle
from astropy.visualization.wcsaxes import SphericalCircle

#seaborn
import seaborn as sns 

#reprojection
from reproject import reproject_interp

#gifs
import imageio 


import pyregion


# -----------------------------------------------------------
# 1. PLOT STYLE
# -----------------------------------------------------------

def styleplots(presentation=False, labelsizes=15, ticksizes=15, txtcolor = 'black', bgcolor = 'white'):
    """
    Make plots pretty and labels clear.
    """
    if presentation==True:
        txtcolor = 'white'
        bgcolor = 'none'
        
    
    plt.rcParams['text.usetex'] = True

    plt.rcParams['font.size'] = labelsizes
    plt.rcParams['axes.titlesize'] = labelsizes
    plt.rcParams['axes.labelsize'] = labelsizes
    plt.rcParams['xtick.labelsize'] = ticksizes
    plt.rcParams['ytick.labelsize'] = ticksizes
    
    #font, font size and color
    # plt.rc('font', family='Helvetica')
    plt.rcParams['font.size'] = labelsizes
    plt.rcParams['text.color'] = txtcolor
    plt.rcParams['axes.labelcolor'] = txtcolor
    plt.rcParams['xtick.color'] = txtcolor
    plt.rcParams['ytick.color'] = txtcolor
    
    #axes 
    plt.rcParams['figure.facecolor'] = bgcolor
    plt.rcParams['axes.facecolor'] = bgcolor
    
    #axes label
    plt.rcParams['axes.titlesize'] = labelsizes
    plt.rcParams['axes.labelsize'] = labelsizes
    plt.rcParams['xtick.labelsize'] = ticksizes
    plt.rcParams['ytick.labelsize'] = ticksizes
    plt.rcParams['axes.labelpad'] = 5
    plt.rcParams['axes.titlepad'] = 5
    
    
    #fig size and saving figures
    
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['figure.figsize'] = [8, 8]
    plt.rcParams['savefig.facecolor'] = bgcolor
    
    return None

# -----------------------------------------------------------
# 2. FIGURE COMPONENTS
# -----------------------------------------------------------

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

def contours(vmax, vmin, nlev):
    nlevels = (vmax - vmin)/nlev
    levels=np.arange(vmin, vmax, nlevels)
    
    return levels


def colorbar(mappable, location="top", pad=0.05, width=5):
    '''
    this includes a hack becuase, if using WCS projection, 
    the astropy WCSAxes object Handles the colorbar ticks 
    through a coordinate map container (you can find it in 
    cbar.ax.coords) instead of the xaxis/yaxis attributes
    See here for more info: https://stackoverflow.com/questions/47060939/matplotlib-colorbar-and-wcs-projection
    '''
    pad = 0.05
    # if location is not "top":
    #     pad = float(input("pad value: enter 0.05 for usual value"))

    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("{}".format(location), size="{}%".format(width), 
                                pad=pad, axes_class=maxes.Axes) # The trick to make sure
    
    if location=="right" or location=="left":
        orientation="vertical"
        cbar = plt.colorbar(mappable, cax=cax, orientation=orientation)
        location = "default"
        cbar.ax.xaxis.set_ticks_position('{}'.format(location))
    
    if location=="top" or location=="bottom":
        #print("here")
        orientation="horizontal"
        cbar = plt.colorbar(mappable, cax=cax, orientation=orientation)
        cbar.ax.xaxis.set_ticks_position('{}'.format(location))
    
    tick_locator = ticker.MaxNLocator(nbins=4)
    cbar.locator = tick_locator
    cbar.update_ticks()
    
    plt.sca(ax)
    return cbar

def pyregion_color(shape, saved_attrs):
    """ selects color for regions drawn """

    attr_list, attr_dict = saved_attrs
    attr_dict["color"] = "green"
    kwargs = pyregion.mpl_helper.properties_func_default(shape, (attr_list, attr_dict))

    return kwargs

def add_beam(ax, hdr, xmin, ymin, xadd = 13, yadd = 10, text=None, 
                bcolor = 'k', fsize=10, fc='none', ytxt_add=0, 
                xtxt_add=0):
    #wcs
    wcs = WCS(hdr, naxis=2)
    
    #beam aperture
    bx, by = xmin + xadd, ymin + yadd #where the beam text will be
    
    if text:
        ax.text(bx+xtxt_add, by+ytxt_add, str(text), ha='center', va='center', 
                color=bcolor, fontsize=fsize)
    
    beam_radec = wcs.wcs_pix2world(bx, by + 6, 0) #where the beam circle will be
    beam_cent = (beam_radec[0], beam_radec[1])
    beam_e = Ellipse(beam_cent, hdr['bmaj'], hdr['bmin'], edgecolor=bcolor, 
                     lw = 2, facecolor=fc, transform=ax.get_transform('fk5'))
    
    ax.add_patch(beam_e)
    return beam_e

def add_circle(ax, xy_center, radius, coord_frame='fk5', lw=2, 
                facecolor='none', ls='--', color='r', alpha=1, 
                label=None, sphere=False):
    """
    This function adds a circle to the given axis
    ax: mpl ax
    xy_center: [xcenter (int), ycenter(int)]
    sphere: set True if making circle on celstial image
    https://docs.astropy.org/en/stable/visualization/wcsaxes/overlays.html
    """
    radius_deg = radius.to(u.deg).value
    xcenter, ycenter = xy_center[0], xy_center[1]

    if sphere is True:
        circle = SphericalCircle((xcenter * u.deg, ycenter * u.deg), radius, 
                            edgecolor=color, lw=lw, facecolor=facecolor, ls=ls, 
                            transform=ax.get_transform(coord_frame), 
                            alpha=alpha, label=label) 

    else:
        circle = plt.Circle((xcenter, ycenter), radius_deg, edgecolor=color, 
                            lw=lw, facecolor=facecolor, ls=ls, 
                            transform=ax.get_transform(coord_frame), 
                            alpha=alpha, label=label) 
    ax.add_patch(circle)

    return None


def hide_plot_labels(ax, hdr, hidex=None, hidey=None, label_only=None):
    """
    This function hides x and y labels. 
    Set hdr = None if not a fits file!
    """
    if label_only:
        if hidex:
            if hdr:
                ax.coords[0].set_axislabel(' ')
            else:
                ax.set_xlabel(' ')

        if hidey:
            if hdr:
                ax.coords[1].set_axislabel(' ')
            else:
                ax.set_ylabel(' ')
        elif hidex is None and hidey is None:
            ax.axis('off')
    else:

        if hidex:
            if hdr:
                ax.coords[0].set_axislabel('')
                ax.coords[0].set_ticklabel_visible(False)
                ax.coords[0].set_ticks_visible(False)
            else:
                ax.set_xlabel('')
                ax.xaxis.set_ticklabels([])
        if hidey:
            if hdr:
                ax.coords[1].set_axislabel('')
                ax.coords[1].set_ticklabel_visible(False)
                ax.coords[1].set_ticks_visible(False)
            else:
                ax.set_ylabel('')
                ax.yaxis.set_ticklabels([])
        
        # elif hidex is None and hidey is None:
        #     ax.axis('off')

    return None

def contour_levels(vmax, vmin, nlev):
    """
    This function creates contour levels
    vmax: float
    vmin: float
    nlev: int
    """
    nlevels = (vmax - vmin)/nlev
    levels=np.arange(vmin, vmax, nlevels)
    
    return levels

def region_shapes(region_string, hdr):
    """
    region_string: string from ds9 when u click list regions
    """

    regions = pyregion.parse(region_string).as_imagecoord(hdr)
    patch_list, artist_list = regions.get_mpl_patches_texts() 
    
    return patch_list

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


def add_scale(ax, scale_dist, hdr, cosmology, scale_angle=None, color='k', lw=1, fsize=12, manual_lim=None, 
                text=None, xtf=0.94, ytf=0.06, ytxt=0.02, hdr_cdelt='CDELT2', zg=0.335):
    #boundaries
    xmin, xmax = ax.get_xlim()[0], ax.get_xlim()[1]
    ymin, ymax = ax.get_ylim()[0], ax.get_ylim()[1]

    if manual_lim:
        xmin, xmax = manual_lim[0][0], manual_lim[0][1]
        ymin, ymax = manual_lim[1][0], manual_lim[1][1]
    
    xpix_dist = (xmax - xmin) * abs(hdr[hdr_cdelt])

    kpc_per_arcsec = cosmology.kpc_proper_per_arcmin(z=zg).to(u.kpc / u.arcsec)

    if scale_dist is not None:
        angle = (scale_dist/kpc_per_arcsec).to(u.arcsec)#.value
        xdistance = (kpc_per_arcsec * xpix_dist * u.deg).to(u.pc)
        scale_factor = (xdistance/scale_dist).cgs
    
    if scale_angle is not None:
        angle = scale_angle
        angle_dist = (kpc_per_arcsec * angle).to(u.kpc)
        scale_factor = ((xpix_dist * u.deg)/angle).cgs


    #scales
    xdelt = (abs(xmax - xmin)/scale_factor).value
    #print(xdelt)

    xs = [xmin + (xmax - xmin)*xtf - xdelt, xmin + (xmax - xmin)*xtf]
    ys = [ymin + (xmax - xmin)*ytf, ymin + (xmax - xmin)*ytf]
    #print(xs, ys)

    ax.plot(xs, ys, color=color, lw=lw)

    if text:
        ax.text(np.mean(xs), ys[0] + (ymax - ys[0]) * ytxt, '{}'.format(text), 
                ha='center', va='center', color=color, fontsize=fsize)

    else:
        if scale_dist is not None:
            ax.text(np.mean(xs), ys[0] + (ymax - ys[0]) * ytxt, '{} ({:.1f}")'.format(scale_dist, angle.value), 
                    ha='center', va='center', color=color, fontsize=fsize)
        if scale_angle is not None:
            ax.text(np.mean(xs), ys[0] + (ymax - ys[0]) * ytxt, '{:.0f} kpc'.format(angle_dist.value), 
                    ha='center', va='center', color=color, fontsize=fsize)
    return ax


# -----------------------------------------------------------
# 3. FITS Handling
# -----------------------------------------------------------

def open_fits(filename, extension=0, naxis=None):
    """ Open FITS File

    Parameters:
        filename (str)
        extension (int): hdu extension
        naxis (int): for WCS

    Return:
        hdu, header, wcs
    """
    
    hdu = fits.open(filename)[extension] # Open the .fits file, return a list of Header/Data Units (HDUs). 
    hdr = hdu.header #The header of the data fits file and it tells you all technical info.
    w = WCS(hdr, naxis=naxis)
    return hdu, hdr, w

def reproject(f1_moment_fits, f2_moment_fits, save_path, name='test.fits'):
    """
    This function opens the HDUs for fits files 1 and 2, picks up their WCSs, 
    drops needless f2 and f1 WCS axes, reprojects f1 to f2, 
    then outputs new HDUs.
    """
    print(name)
    f1HDU = fits.open(f1_moment_fits)
    f2HDU = fits.open(f2_moment_fits)

    w_f2 = WCS(f2HDU[0])
    w_f1 = WCS(f1HDU[0])
    
    f1_shape = w_f1.pixel_shape
    f2_shape = w_f2.pixel_shape
    
    
    if w_f1.pixel_n_dim == 3:
        f1_2wcs = w_f1.dropaxis(2)
    elif w_f1.pixel_n_dim == 2:
        f1_2wcs = w_f1
        
    if w_f2.pixel_n_dim == 3:
        f2_2wcs = w_f2.dropaxis(2)
    elif w_f2.pixel_n_dim == 2:
        f2_2wcs = w_f2

    newf2Header = f2_2wcs.to_header()
    newf2HDU = fits.PrimaryHDU(data=f2HDU[0].data, header=newf2Header)
    
    newf1Header = f1_2wcs.to_header()
    newf1HDU = fits.PrimaryHDU(data=f1HDU[0].data, header=newf1Header)
    
    #reproject f1 onto f2
    
    registered_f1_data, registered_f1_footprint = reproject_interp(newf1HDU, newf2Header, 
                                                                    shape_out=(f2_shape[1],
                                                                     f2_shape[0]))

    fig = plt.figure(figsize=(12,8))
    
    ax1 = plt.subplot(1,2,1, projection=WCS(newf1Header))
    ax1.imshow(f1HDU[0].data, origin='lower', norm=LogNorm())
    ax1.coords['ra'].set_axislabel('RA')
    ax1.coords['dec'].set_axislabel('Dec')
    ax1.set_title('F1 Original')

    ax2 = plt.subplot(1,2,2, projection=WCS(newf2Header))
    ax2.imshow(registered_f1_data, origin='lower')
    ax2.coords['ra'].set_axislabel('RA')
    ax2.coords['dec'].set_axislabel('Dec')
    ax2.set_title('Registered F1 Image')
    
    plt.show()
    print("Successfully reprojected {} to {}".format(f2_moment_fits, f1_moment_fits))
    
    header = w_f2.to_header()
    
    print("saving new header to products folder")
    
    
    filename = ''.join((save_path, 'reprojected_{}'.format(name)))

    fits.writeto(filename, registered_f1_data, header, overwrite=True)
    
    return registered_f1_data, header, filename

def zoom_fits(hdr, xcoord, ycoord, radius, spitzer=False, naxis=2):
    """
    This function zooms in on a certain area of a fits map
    radius: with astropy units
    xcoord: ra/glon
    ycoord:dec/glat
    """
    w = WCS(hdr, naxis=naxis)

    rad_deg = radius.to(u.deg).value
    if spitzer is True:
        pixel_size = hdr['CD2_2']
    else:
        pixel_size = hdr['CDELT2']
    pixel_radius = int(rad_deg/pixel_size)
    #print(pixel_radius, rad_deg, pixel_size)

    if w.naxis==2:
        pixel_coords = w.wcs_world2pix(xcoord, ycoord, 1) #coordinate location in pixels
        px, py = int(pixel_coords[1]), int(pixel_coords[0]) 
    
    if w.naxis==3:
        pixel_coords = w.wcs_world2pix(xcoord, ycoord, 0, 1) #coordinate location in pixels
        px, py = int(pixel_coords[1]), int(pixel_coords[0]) 
    
    if w.naxis==4:
        pixel_coords = w.wcs_world2pix(xcoord, ycoord, 0, 0, 1) #coordinate location in pixels
        px, py = int(pixel_coords[1]), int(pixel_coords[0]) 

    xmin, xmax = px - pixel_radius, px + pixel_radius
    ymin, ymax = py - pixel_radius, py + pixel_radius
    xrange, yrange = [xmin, xmax], [ymin, ymax]

    return xrange, yrange

def vel_to_pixel(hdu, vel, opp=False, vp=None):
    """
    This function converts velocity in (km/s) to the corresponding
    slice in pixel space
    """
    hdr = hdu.header
    w = WCS(hdr)

    if opp is False:
        #print("not here")
        pixel_coords = w.wcs_world2pix(0, 0, vel * 1000, 1) 
        vel_pixel = int(pixel_coords[2]) - 1
        #print(vel_pixel)

    if opp is True:
        #print("here")
        world_coords = w.wcs_pix2world(0, 0, vp, 1) 
        vel_pixel = world_coords[2]/1e3 
       #print(vel_pixel)

    return vel_pixel

# -----------------------------------------------------------
# 4. PLOT FITS MAPS
# -----------------------------------------------------------

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
    # ax.grid(False)
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


def ax_plot(ax, data, hdr, cmap, clip=None, contour=False, contour_levels=10, vmin=None, 
            vmax=None, units=None, bcolor='white', xy_coords = None, abeam=False, 
            alpha = 1, cb_color='k', location="top", cbar_no_ticks=False, 
            cb_width=8, cb_pad=-38, cb_fsize=12, cb_txt_loc='center', show_cbar=True,
            cb_txt=None, cb_txt_pos=[0, 0], contour_lines=False):
    
    if clip is not None: 
        nanmask = (data < clip)
        data[nanmask] = np.nan
    
    #plot image
    if contour is False:
        im = ax.imshow(data, cmap=cmap, vmin = vmin, vmax = vmax, alpha=alpha)
    
    elif contour is True:
        levels = contours(vmax, vmin, contour_levels)
        im = ax.contourf(data, levels=levels, cmap=cmap, extend='min', alpha = alpha)

        for c in im.collections:
            if contour_lines == False:
                c.set_edgecolor("face")
                c.set_rasterized(True)

    # #labels
    # labels(ax, hdr)
    
    #color bar
    cbar = colorbar(im, location=location, width=cb_width)
    cbar.ax.tick_params(direction='in', length=8)
    if units:
        cbar.set_label('%s' % (units), color=cb_color, labelpad=cb_pad, 
                        loc=cb_txt_loc, fontsize=cb_fsize)
    if cb_txt:
        cbar.ax.text(cb_txt_pos[0], cb_txt_pos[1], cb_txt, 
                    color=cb_color, fontsize=cb_fsize)
    if location is 'left' or location is 'right':
        cbar.ax.get_yticklabels()[0].set_visible(False)
        cbar.ax.get_yticklines()[0].set_visible(False)
        cbar.ax.get_yticklabels()[-1].set_visible(False)
        cbar.ax.get_yticklines()[-1].set_visible(False)

    if location is 'top' or location is 'bottom':
        cbar.ax.get_xticklabels()[0].set_visible(False)
        cbar.ax.get_xticklines()[0].set_visible(False)
        cbar.ax.get_xticklabels()[-1].set_visible(False)
        cbar.ax.get_xticklines()[-1].set_visible(False)

    if cbar_no_ticks is True:
        cbar.set_ticks([])
    if show_cbar is False:
        cbar.remove()
    #zoom in
    if xy_coords:
        xmin, xmax, ymin, ymax = xy_coords[0][0], xy_coords[0][1], xy_coords[1][0], xy_coords[1][1]
        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)

        #beam aperture
        if abeam is True:
            beam = add_beam(ax, hdr, xmin, ymin, text='Beam', bcolor = bcolor)
    
    return None

def plot_fits_map(data, hdr, cmap, wcs, filename=None, vmin=None, vmax=None, 
                    figsize=(10,8), slice=0, gmos=None, contour=False, 
                    contour_levels=10, clip=None, clip_max=None, location="top",
                    width=5):    
    """Plot One FITS Map"""
    if hdr:
        wcs = WCS(hdr)
        if gmos:
            wcs = wcs.dropaxis(2)
        
    if filename:
        data, hdr, wcs = open_fits(filename)
        
        
    fig = plt.figure(figsize=figsize)
    if wcs.naxis==3:
        ax = fig.add_subplot(111, projection=wcs, slices=('x', 'y', slice))
        data=data[slice, :, :]
    else:
        ax = fig.add_subplot(111, projection=wcs)

    if clip is not None: 
        nanmask = (data < clip)
        data[nanmask] = np.nan

    if clip_max is not None: 
        nanmask = (data > clip)
        data[nanmask] = np.nan
    
    if contour is False:
        im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax)
    
    elif contour is True:
        levels = np.linspace(vmin, vmax, contour_levels)
        im = plt.contourf(data, levels=levels, cmap=cmap, extend='min')
    
    cbar = colorbar(im, location=location, width=width)
    cbar.ax.tick_params(direction='in', length=8)
    
    if hdr:
        labels(ax, hdr)
    
    return fig, ax

def plot_all_moments(moment_filenames, header, ra=232.7936938, dec=34.2404172, radius=2.5*u.arcsec, 
                    vmin=0, vmax=None, cmap=sns.color_palette("rocket_r", as_cmap=True), clip=None, 
                    units=None, bcolor='k', contour=False, nlevels=10, region_string=None, ncols=5, 
                    nrows=1):
    """ 
    Plot all moments(0,1,2,8,9) for each file given

    moment_filenames: list of filenames
    header: header for each file in moment_filenames
    """
    w = WCS(header)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5*ncols, 5*nrows),
                                subplot_kw={'projection':w, 'slices':('x', 'y')})
    if ncols==1:
        axs = [axes]
    else:
        axs = axes.ravel()

    for idf, filename in enumerate(moment_filenames):

        hdu, header, wcs = open_fits(filename, naxis=2)
        zoom_coords=[ra, dec, radius*3]

        ax = axs[idf]

        #plot moment 
        ax, cb = ax_fits_map(hdu, vmin=vmin, vmax=vmax, cmap=cmap, 
                     coords=zoom_coords, ax=ax, contour=contour, 
                     nlevels=nlevels, clip=clip)
        add_beam(ax, header, 109, 87, bcolor=bcolor, text='Beam')
        labels(ax, header)

        if region_string is not None:
            patch_list = region_shapes(region_string, header)
            for p in patch_list:
                ax.add_patch(p)
    
    return fig

def ax_fits_map(hdu, cmap, ax, filename=None, data=None, hdr=None, factor=None, clip=None, lw=1, vel=0, 
                vslice=None, vmin=None, vmax=None, drop_axis=None, contour=False, nlevels=10, coords=None, 
                spitzer=False, location="top", cbar=True, pad=0.05, clip_max=None, wcs=None, 
                cbar_no_ticks=False, cb_width=5, cb_pad=-38, cb_fsize=12, units=None, cb_color='k',
                plot_contours=False, pc_vmin=0, pc_vmax=1, pc_nlev=10, pc_color='k', pc_lw=1, 
                pc_alpha=1, pc_ls='-', naxis=2, norm=False):   

    """
    This function plots a singular fits map on a given axis
    hdu: fits hdu
    cmap: selected colormap
    wcs: fits wcs info
    ax: figure axis 
    coords: center on this area
    
    give filename if file not too big and set other required variables to none

    drop_axis: wcs axis to drop (int) 
    coords: [ra, dec, radius]  or [glon, glat, radius]
    """
    if hdu:
        hdr = hdu.header
        data = hdu.data
        wcs = WCS(hdr)
        if drop_axis:
            wcs = wcs.dropaxis(drop_axis)
        
    if filename:
        hdu, hdr, wcs = open_fits(filename)
        data = hdu.data #assuming 2d

    if data is not None:
        if wcs is None:
            wcs=WCS(hdr)

    if wcs.naxis==3:
        if vslice is None:
            vslice = vel_to_pixel(hdu, vel)
            if vslice<0:
                vslice=0
            if vslice > hdr['NAXIS3']:
                vslice=hdr['NAXIS3'] - 1
        data = data[vslice, :, :]
    
    if wcs.naxis==4:
        data = data[0,0, :, :]

    if factor:
        data = data/factor

    if clip is not None:
        data = data.astype('float32')
        nanmask = (data < clip)
        data[nanmask] = np.nan

    if clip_max is not None: 
        data[np.where(data > clip_max)] = np.nan
    
    if contour is not True:
        im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    
    elif contour is True:
        
        levels = contour_levels(vmax, vmin, nlevels)
        im = plt.contourf(data, levels=levels, cmap=cmap, extend='min', aspect='auto')
        for c in im.collections:
            c.set_edgecolor("face")
            c.set_rasterized(True)
    
    if norm is True:
        im.set_norm(LogNorm(vmin=vmin, vmax=vmax))
    if plot_contours is True:
        overlay(ax,  hdu=None, data=data, hdr=hdr, alpha=pc_alpha, vmin=pc_vmin, vmax=pc_vmax, nlev=pc_nlev, contours=True, 
                cont_color=pc_color, lw=pc_lw, ls=pc_ls)
    
    if hdr:
        #print(ax)
        labels(ax, hdr)
        
    if coords:
        ycoords, xcoords = zoom_fits(hdr, coords[0], coords[1], coords[2], spitzer=spitzer, naxis=naxis)

        ax.set_xlim(xcoords[0], xcoords[1])
        ax.set_ylim(ycoords[0], ycoords[1])
        #print("here, x:{}, y: {}".format(xcoords[0], ycoords[0]))
    
    cb = colorbar(im, location=location, pad=pad, width=cb_width)
    cb.ax.tick_params(direction='in', length=8, labelsize=cb_fsize)

    if units:
        cb.set_label('%s' % (units), color=cb_color, labelpad=cb_pad, 
                        loc='center', fontsize=cb_fsize)

    if location is 'left' or location is 'right':
        cb.ax.get_yticklabels()[0].set_visible(False)
        cb.ax.get_yticklines()[0].set_visible(False)
        cb.ax.get_yticklabels()[-1].set_visible(False)
        cb.ax.get_yticklines()[-1].set_visible(False)

    if location is 'top' or location is 'bottom':
        cb.ax.get_xticklabels()[0].set_visible(False)
        cb.ax.get_xticklines()[0].set_visible(False)
        cb.ax.get_xticklabels()[-1].set_visible(False)
        cb.ax.get_xticklines()[-1].set_visible(False)

    if cbar_no_ticks is True:
        cb.set_ticks([])

    if cbar is False:
        cb.remove()

    cb.solids.set_rasterized(True)
    cb.solids.set_edgecolor("face")

    return ax, cb


# -----------------------------------------------------------
# 5. MOVIE
# -----------------------------------------------------------

def alma_movie(cube_file, num_chans,vmin, vmax, movie_dir, fps, ap_type=[0], 
                xy_coords = [110, 205, 107, 185], bcolor='white', hst_data=None, 
                cvmax=None, cmap=cm.magma):
    """
    #num_chans: how long do we want each movie
    #vmin,vmax: velocity range of velocity channels
    #movie_dir: where to place the finished gifs
    #fps: how many frames per second. 2 is a good value
    #ap_type: list of apertures to be displayed? [0] = total, [1] = north, [2] = south or [0,1] etc
    hst_data: loaded hdu
    """
    
    #make directories for final gifs and gif images
    os.makedirs(movie_dir, exist_ok=True)
    img_dir = movie_dir + 'gif_images/'
    os.makedirs(img_dir, exist_ok=True)
    
    #remove any previous images
    for f in glob.glob(img_dir + "*.png"):
        os.remove(f)

    img_files = []
    
    #system info
    ra = 232.7936938
    dec = 34.2404172
    radius = Angle(2.5, u.arcsec)
    rad_deg = radius.to(u.deg)
    
    #read in data
    hdu = fits.open(cube_file)[0]
    hdr = hdu.header
    w = WCS(hdr)
    
    pvmin, pvmax = int(w.wcs_world2pix(ra, dec, vmin*1000, 1)[2]), int(w.wcs_world2pix(ra, dec, vmax*1000, 1)[2])
    slices = np.linspace(pvmin, pvmax, num_chans)
    slices = [ int(x) for x in slices]
        
    png_files = []
    
    #plot aesthetics
    cmap = cm.magma
        
    for idx, slicex in enumerate(slices):

        #plotting
        fig, ax = plot_fits_map(hdu.data, hdu.header, cmap=cmap, vmin=0.0005, vmax=cvmax, 
                                wcs=None, slice=slicex) 
        
        xmin, xmax, ymin, ymax = 110, 205, 107, 185
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

        #add apertures 
        rad = rad_deg.value
        names = ['Total', 'North', 'South']
        ras = [ra+10e-5,  ra-1e-5, ra+2e-4]#+2.5e-5 #ra
        decs = [dec-3.8e-4, dec-6e-5, dec-7e-4] #dec
        wds = [rad*2.4, rad, rad*1.2] #width
        hts = [rad*1.2, rad * 1.1, rad*1.1] #height
        angs = [290, 10 , 0 ]#angle
        colors = ['yellow', 'lightgrey', 'red'] #edgecolor

        for ide in ap_type:
            era, edec, ewd, eht, eang, ecolor, ename = ras[ide], decs[ide], wds[ide], hts[ide], angs[ide], colors[ide], names[ide]
            aperture = Ellipse((era, edec), ewd, eht, angle = eang,  edgecolor=ecolor, facecolor='none', 
                       lw=1.2, ls = '--', transform=ax.get_transform('fk5'), label = '%s' % (ename))

            ax.add_patch(aperture)
            
        #add beam aperture
        add_beam(ax, hdr, xmin, ymin, text='Beam', bcolor = 'w')
        
        
        #add labels
        labels(ax, hdr)
        
        #hst contours
        if hst_data:
            #f606 = fits.open(hst_file)[0]
            overlay(ax, hst_data, 1, vmin=0.035, vmax=0.3, contours=True, cont_color=bcolor)


        #ax.set_title('ALMA CO(3-2)a \n')
        
        xlim = int((ax.get_xbound()[1])/1.2)
        ylim = int((ax.get_ybound()[1]) - (0.05*ax.get_ybound()[1]))
        vel = w.wcs_pix2world(0,0,slicex, 0)[2]/1000
        ax.text(xlim, ylim, '%s km/s' % (round(vel, 1)), style='italic', 
                bbox={'facecolor': 'white', 'alpha': 0.9, 'pad': 10}, color = 'black')

        #save images
        img_name = img_dir + '{}'.format(idx) + '.png'
        #plt.rcParams['savefig.facecolor'] = None
        fig.savefig(img_name, bbox_inches='tight', dpi=100, pad_inches=0)
        img_files.append(img_name)
        

        plt.close(fig) # don't spam me with a gajillion figures

    #make the gif
    movie_name = 'SDSS_J1531+3414_{}_movie'.format(ap_type)
    gif_name = movie_dir + movie_name + '.gif'
    gif_frames = []

    # Remove any old GIFs you might have made
    if os.path.isfile(gif_name):
        os.remove(gif_name)

    for filename in img_files:
        gif_frames.append(imageio.imread(filename))

    imageio.mimsave(gif_name, gif_frames, fps=fps)
    
    return None
