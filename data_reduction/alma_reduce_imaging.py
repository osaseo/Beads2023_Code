

########################################
# Check CASA version

import re
import os

if casadef.casa_version < '4.4.0' :
    sys.exit("Please use CASA version greater than or equal to 4.4.0 with this script")


##################################################
# Create an Averaged Continuum MS


finalvis='calibrated_final.ms' # This is your output ms from the data
                               # preparation script.

# Use plotms to identify line and continuum spectral windows.
if not os.path.exists('calibrated_final_cont.ms'):
	plotms(vis=finalvis, xaxis='channel', yaxis='amplitude',
	       ydatacolumn='data',
	       avgtime='1e8', avgscan=True, avgchannel='1', 
	       iteraxis='spw' )





	# Set spws to be used to form continuum
	contspws = '0,1,2,4,5,6,8,9,10'

	# If you have complex line emission and no dedicated continuum
	# windows, you will need to flag the line channels prior to averaging.
	flagmanager(vis=finalvis,mode='save',
	            versionname='before_cont_flags')

	initweights(vis=finalvis,wtmode='weight',dowtsp=True)

	# Flag the "line channels"
	flagchannels='2:25~40;100~127,6:25~40;100~127,10:25~40;100~127' 

	flagdata(vis=finalvis,mode='manual',
	          spw=flagchannels,flagbackup=False)

	# check that flags are as expected, NOTE must check reload on plotms
	# gui if its still open.

	plotms(vis=finalvis,yaxis='amp',xaxis='channel',
	       avgchannel='1',avgtime='1e8',avgscan=True,iteraxis='spw') 

	# Average the channels within spws
	contvis='calibrated_final_cont.ms'
	rmtables(contvis)
	os.system('rm -rf ' + contvis + '.flagversions')


	split2(vis=finalvis,
	     spw=contspws,      
	     outputvis=contvis,
	       width=[8,8,8,8,8,8,8,8,8], 
	     datacolumn='data')


	# Check the weights. You will need to change antenna and field to
	# appropriate values
	plotms(vis=contvis, yaxis='wtsp',xaxis='freq',spw='',antenna='DA42',field='5')

	# If you flagged any line channels, restore the previous flags
	flagmanager(vis=finalvis,mode='restore',
	            versionname='before_cont_flags')

	# Inspect continuum for any problems
	plotms(vis=contvis,xaxis='uvdist',yaxis='amp',coloraxis='spw')

# #############################################
# Image Parameters


# source parameters
# ------------------

field='5' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
imagermode='csclean' # uncomment if single field 
# imagermode='mosaic' # uncomment if mosaic or if combining one 7m and one 12m pointing.
# phasecenter=3 # uncomment and set to field number for phase
                # center. Note lack of ''.  Use the weblog to
                # determine which pointing to use. Remember that the
                # field ids for each pointing will be re-numbered
                # after your initial split. You can also specify the
                # phase center using coordinates, e.g.,
                # phasecenter='J2000 19h30m00 -40d00m00'

# image parameters.
# ----------------

cell='0.1arcsec' # cell size for imaging.
imsize = [300,300] # size of image in pixels.

# velocity parameters
# -------------------

outframe='lsrk' # velocity reference frame. See science goals.
veltype='radio' # velocity type. See note below.


# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean. 

weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.1mJy'

#############################################
# Imaging the Continuuum

# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'         
contimagename = 'calibrated_final_cont'

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=contvis)
#delmod(vis=contvis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=contvis,
      imagename=contimagename,
      field=field,
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      interactive = False,
      imagermode = imagermode)


# rms is ~16 uJy

# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an optional step.
#contmaskname = 'cont.mask'
##rmtables(contmaskname) # if you want to delete the old mask
#os.system('cp -ir ' + contimagename + '.mask ' + contmaskname)

##############################################
# Image line emission


finalvis = 'calibrated_final.ms'
linevis = finalvis # uncomment if you neither continuum subtracted nor self-calibrated your data.
# linevis = finalvis + '.contsub' # uncomment if continuum subtracted
# linevis = finalvis + '.contsub.selfcal' # uncommment if both continuum subtracted and self-calibrated
# linevis = finalvis + '.selfcal' # uncomment if just self-calibrated (no continuum subtraction)

sourcename ='SDSS_J15' # name of source
linename = 'CO32_40kms' # name of transition (see science goals in OT for name) 
lineimagename = sourcename+'_'+linename # name of line image

# CO(3-2) rest frequency is 345.796 GHz. z = 0.335
# so 345.796 / (1 + 0.335) = 259.02322 GHz

restfreq='259.02322GHz'
# THIS WAS WHAT WHAS ORIGINALLY THERE
#restfreq='259.04269GHz' # Typically the rest frequency of the line of
                        # interest. If the source has a significant
                        # redshift (z>0.2), use the observed sky
                        # frequency (nu_rest/(1+z)) instead of the
                        # rest frequency of the
                        # line.

spw='3,7,11' # uncomment and replace with appropriate spw if necessary.
threshold = '0.45mJy'

start='-600km/s' # start velocity. See science goals for appropriate value.
width='40km/s' # velocity width. See science goals.
nchan = 40  # number of channels. See science goals for appropriate value.


# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=linevis)
#delmod(vis=linevis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(lineimagename + ext)

clean(vis=linevis,
      imagename=lineimagename, 
      field=field,
      spw=spw,
      mode='velocity',
      start=start,
      width=width,
      nchan=nchan, 
      outframe=outframe, 
      veltype=veltype, 
      restfreq=restfreq, 
      niter=niter,  
      threshold=threshold, 
      interactive=False,
      cell=cell,
      imsize=imsize, 
      weighting=weighting, 
      robust=robust,
      imagermode=imagermode)


# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an
# optional step.
# linemaskname = 'line.mask'
## rmtables(linemaskname) # uncomment if you want to overwrite the mask.
# os.system('cp -ir ' + lineimagename + '.mask ' + linemaskname)

##############################################
# Apply a primary beam correction

import glob

myimages = glob.glob("*.image")

rmtables('*.pbcor')
for image in myimages:
    pbimage = image.rsplit('.',1)[0]+'.flux'
    outfile = image.rsplit('.',1)[0]+'.pbcor'
    impbcor(imagename=image, pbimage=pbimage, outfile = outfile)

##############################################
# Export the images

import glob

myimages = glob.glob("*.pbcor")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*.flux")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True) 

##############################################
# Create Diagnostic PNGs


#os.system("rm -rf *.png")
#mycontimages = glob.glob("calibrated*.image")
#for cimage in mycontimages:
#    max=imstat(cimage)['max'][0]
#    min=-0.1*max
#    outimage = cimage+'.png'
#    os.system('rm -rf '+outimage)
#    imview(raster={'file':cimage,'range':[min,max]},out=outimage)


# this will have to be run for each sourcename
#sourcename='' # insert source here, if it isn't already set
#mylineimages = glob.glob(sourcename+"*.image")
#for limage in mylineimages:
#    rms=imstat(limage,chans='1')['rms'][0]
#    mom8=limage+'.mom8'
#    os.system("rm -rf "+mom8)
#    immoments(limage,moments=[8],outfile=mom8)
#    max=imstat(mom8)['max'][0]
#    min=-0.1*max
#    os.system("rm "+mom8+".png")
#    imview(raster={'file':mom8,'range':[min,max]},out=mom8+'.png')


##############################################
# Analysis

# For examples of how to get started analyzing your data, see
#     https://casaguides.nrao.edu/index.php/TWHydraBand7_Imaging_4.3
#     
