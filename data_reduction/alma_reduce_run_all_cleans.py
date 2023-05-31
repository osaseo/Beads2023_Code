# Script must be run in CASA >6.2.1.7 as execfile('run_all_cleans.py')

import time
import numpy as np

import matplotlib.pyplot as plt
from astropy.io import fits

# Check CASA version
try:
    import casalith
except:
    print("Script requires CASA 6.0 or greater")

if casalith.compare_version("<", [6, 2, 1, 7]):
    print("Please use CASA version greater than or equal to 6.2.1.7 with this script")


def calc_num_channels(start=-800, width=40):
    ''' Let's just use symmetric start/end points around the systemic
    '''
    nchan = int((abs(start) * 2) / width)
    return nchan


def save_fits(run_output_name):
    exportfits(imagename=run_output_name+'.image',
               fitsimage=run_output_name+'.fits', overwrite=True)


def plot_results(run_output_name):

    filename_extensions = ['.image', '.mask', '.model',
                           '.pb', '.psf', '.residual']

    ff, aa = plt.subplots(2, 3, figsize=(18, 12))
    for ii, extension in enumerate(filename_extensions):
        exportfits(imagename=run_output_name+extension,
                   fitsimage=run_output_name+extension+'.fits', overwrite=True)
        xx, yy = int(ii/3), ii % 3
        im = aa[xx, yy].imshow(fits.getdata(
            run_output_name+extension+'.fits')[0, 0, :, :])
        plt.colorbar(im, ax=aa[xx, yy])
        aa[xx, yy].set_title(extension)
    plt.suptitle('{}'.format(run_output_name), fontsize=18)

    plt.savefig(run_output_name+'_quicklook_'+'.png')
    plt.close()


def run_tclean(measurement_sets, run_output_name, width, weighting, uvtaper):

    tclean(vis=measurement_sets,
           imagename=run_output_name,
           field='5',
           spw='3,7,11',
           specmode='cube',  # appropriate for non-ephemeris source
           perchanweightdensity=True,
           start='-700km/s',
           width=str(width)+'km/s',
           nchan=calc_num_channels(width=width),
           outframe='lsrk',
           veltype='optical',
           restfreq='259.02322GHz',  # so 345.796 / (1 + 0.335) = 259.02322 GHz
           niter=5000,
           threshold='0.45mJy',
           interactive=False,
           cell='0.1arcsec',
           imsize=[300, 300],
           weighting=weighting,
           robust=0.5,
           uvtaper=uvtaper,
           gridder='standard',  # for a single field
           pbcor=True,
           restoringbeam='common',
           usepointing=False,
           parallel=True)

    print('Done with {}, exporting to FITs...'.format(run_output_name))
    exportfits(imagename=run_output_name+'.image.pbcor', fitsimage=run_output_name +
               '.image.pbcor.fits', velocity=True, optical=True, overwrite=True)

    print("Exported {}".format(run_output_name+'.image.pbcor.fits'))

    plot_results(run_output_name)


startTime = time.time()

measurement_sets = ['647_calibrated_final.ms', '649_calibrated_final.ms']
weightings_to_use = ['natural', 'briggs']
widths_to_use = [10, 20, 40, 80]
tapers_to_use = ['0.6arcsec', '0.8arcsec', '1.0arcsec', '1.2arcsec']

for weight in weightings_to_use:
    print('Running {} run...'.format(weight))
    for width in widths_to_use:
        print('Running {}km/s NOTAPER run...'.format(width))
        run_name = 'SDSS1531_CO32_{}kms_{}_notaper'.format(width, weight)
        run_tclean(measurement_sets, run_output_name=run_name,
                   width=width, weighting=weight, uvtaper=[])

        for taper in tapers_to_use:
            print('Running {}km/s TAPER run with {} taper...'.format(width,taper))
            run_name = 'SDSS1531_CO32_{}kms_{}_taper_{}'.format(
                width, weight, taper)
            run_tclean(measurement_sets, run_output_name=run_name,
                       width=width, weighting=weight, uvtaper=taper)

execution_time = np.round((time.time() - startTime) / 3600, 2)
print('Finished in {} hours. Phew.'.format(execution_time))
