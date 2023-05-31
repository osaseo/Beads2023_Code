# run this code in an active CIAO environment in a directory WITHOUT spaces in the path

#download Chandra data
download_chandra_obsid 17218,18689

#reprocess Chandra data
punlearn ardlib
chandra_repro 17218 verbose=5
punlearn ardlib
chandra_repro 18689 verbose=5


# merge_obs simply combines the reprojet_obs and flux_obs scripts
#merges obsids 17218 and 18689
punlearn ardlib
punlearn merge_obs
merge_obs "*/repro/*evt*[ccd_id=7]" merged/beads_xray_bin0.5 bin=0.5 bands=broad,csc verbose=5
merge_obs "*/repro/*evt*[ccd_id=7]" merged/beads_xray_bin1 bin=1 bands=broad,csc verbose=5
merge_obs "*/repro/*evt*[ccd_id=7]" merged/beads_xray_bin2 bin=2 bands=broad,csc verbose=5
merge_obs "*/repro/*evt*[ccd_id=7]" merged/beads_xray_bin4 bin=4 bands=broad,csc verbose=5


