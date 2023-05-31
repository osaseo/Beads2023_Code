#Run this script in an active CIAO environment 

#To find evidence of cavities, we create unsharp masks in CIAO for a variety of different bins and smoothed with different gaussians

# We create masks by dividing the image smoothed with a certain radius Gaussian by one smoothed with another radius Gaussian. 
# The radii are given below:
# - 0.98" - 9.8"
# - 1.5" - 8"
# - 1.5" - 10"
# - 2" - 10"
# - 0.7" - 6"
# - 1" - 8"
# - 1" - 6"

#This script only creates them for the bin1 images, but you can modify this script to create them for the other bins if you're curious!

#create the gaussian smoothed images
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g098.img "lib:gaus(2,5,1,0.98,0.98)" meth=fft mode=h clob+
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g98.img "lib:gaus(2,5,1,9.8,9.8)" meth=fft mode=h clob+

#create unsharp mask
dmimgcalc bin1_broad_flux_g098.img bin1_broad_flux_g98.img bin1_broad_flux_unsharpmask_098_98.img sub clob+

#for the rest

#1.5" - 8"
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g15.img "lib:gaus(2,5,1,1.5,1.5)" meth=fft mode=h clob+
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g8.img "lib:gaus(2,5,1,8.0,8.0)" meth=fft mode=h clob+
dmimgcalc bin1_broad_flux_g15.img bin1_broad_flux_g8.img bin1_broad_flux_unsharpmask_15_8.img sub clob+

#1.5" - 10"
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g10.img "lib:gaus(2,5,1,10.0,10.0)" meth=fft mode=h clob+
dmimgcalc bin1_broad_flux_g15.img bin1_broad_flux_g10.img bin1_broad_flux_unsharpmask_15_10.img sub clob+

#2" - 10"
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g2.img "lib:gaus(2,5,1,2.0,2.0)" meth=fft mode=h clob+
dmimgcalc bin1_broad_flux_g2.img bin1_broad_flux_g10.img bin1_broad_flux_unsharpmask_2_10.img sub clob+

#0.7" - 6"
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g07.img "lib:gaus(2,5,1,0.7,0.7)" meth=fft mode=h clob+
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g6.img "lib:gaus(2,5,1,6.0,6.0)" meth=fft mode=h clob+
dmimgcalc bin1_broad_flux_g07.img bin1_broad_flux_g6.img bin1_broad_flux_unsharpmask_07_6.img sub clob+

#1" - 8"
aconvolve beads_xray_bin1_broad_flux.img bin1_broad_flux_g1.img "lib:gaus(2,5,1,1.0,1.0)" meth=fft mode=h clob+
dmimgcalc bin1_broad_flux_g1.img bin1_broad_flux_g8.img bin1_broad_flux_unsharpmask_1_8.img sub clob+

#1" - 6"
dmimgcalc bin1_broad_flux_g1.img bin1_broad_flux_g6.img bin1_broad_flux_unsharpmask_1_6.img sub clob+