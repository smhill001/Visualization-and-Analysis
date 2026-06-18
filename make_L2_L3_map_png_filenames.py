def make_L2_L3_map_png_filenames(obskey,filename,Level,LonSys,LatLims,LonLimsEast,
                                 coef=0.0,FiveMicron=False,param_name=False,
                                 dataversion=2):
    import numpy as np
    import make_lat_lon_str as MLLS

    if coef==0.0:
        correction='_C0'
    else:
        correction='_C1'

    latstr,lonstr=MLLS.make_lat_lon_str(LatLims,LonLimsEast)
    if FiveMicron:
        fnskeleton=correction+'_Sys'+LonSys+'_N'+\
                    str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                    '_Lon'+str(np.mod(360-LonLimsEast[1],360)).zfill(3)+'-'+\
                        str(np.mod(360-LonLimsEast[0],360)).zfill(3)+'_5micron.png'
    else:
        fnskeleton=correction+'_Sys'+LonSys+'_'+latstr+'_'+lonstr+'.png'

    ###!!!! the file extension and L3 header FILENAME does!
    if Level=='L2':
        fnout=filename+fnskeleton
    elif Level=='L3':
        fnout=obskey+'_SCT_Jupiter'+fnskeleton
        fnout=fnout.replace('.png',' '+param_name+'.png')
    if dataversion=='H':
        fnout=obskey+'_HST_Jupiter'+fnskeleton
        fnout=fnout.replace('.png',' '+param_name+'.png')
    return fnout