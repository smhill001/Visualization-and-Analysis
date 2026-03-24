from config_VA import config_VA

import numpy as np
from scipy.signal import convolve2d
import plot_patch as pp
import matplotlib.pyplot as pl

def rolling_corr2d(A, B, wdeg,dataversion=2):
    
    if dataversion=='H':
        w=int(wdeg*20)
    else:
        w=int(wdeg*2)

    kernel = np.ones((w, w))
    N = w * w

    sumA  = convolve2d(A, kernel, mode='same', boundary='symm')
    sumB  = convolve2d(B, kernel, mode='same', boundary='symm')

    sumA2 = convolve2d(A*A, kernel, mode='same', boundary='symm')
    sumB2 = convolve2d(B*B, kernel, mode='same', boundary='symm')

    sumAB = convolve2d(A*B, kernel, mode='same', boundary='symm')

    numerator = N*sumAB - sumA*sumB

    denomA = N*sumA2 - sumA**2
    denomB = N*sumB2 - sumB**2

    denom = np.sqrt(denomA * denomB)

    r = numerator / denom

    return r

def residual_2d(A,B):
    
    mean_A = np.mean(A)
    amp_A = np.max(np.abs(A - mean_A))
    norm_A = (A - mean_A) / amp_A
    
    mean_B = np.mean(B)
    amp_B = np.max(np.abs(B - mean_B))
    norm_B = (B - mean_B) / amp_B
    
    resid_AB=norm_A-norm_B
    
    return norm_A,norm_B,resid_AB


def L3_Jup_Map_Plot_V2(obskey="20251016UTa",target="Jupiter",
                        CoLatLims=[45,135],LonRng=45,
                        CMpref='subobs',LonSys='2',showbands=False,
                        coef=[0.,0.],subproj='',figxy=[8.0,4.0],FiveMicron=False,
                        plotoptions=["contours","surface"],
                        ROI=False,
                        LimbCorrection=False,dataversion=2,smoothcont=0):
    """
    Created on Sun Nov  6 16:47:21 2022
    
    PURPOSE: Create MAPS of environmental parameters paired with RGB context
             maps. Based on Retrieve_Jup_Atm_2022_P3, which ALSO performed
             the calibration phase. So now I've separated that module into 
             a calibration module, make_L3_env_data.py and this plotting
             module.
             
    EXAMPLES:
        Map_Jup_Atm_P3(obskey="20240925UTa",
                                CoLatLims=[45,135],LonRng=45,
                                CMpref='subobs',LonSys='2',showbands=False,
                                coef=[0.,0.],subproj='',figxy=[8.0,4.0],
                                FiveMicron=False)
        
        Map_Jup_Atm_P3(obskey="20240730UTa",
                                CoLatLims=[45,135],LonRng=45,
                                CMpref='subobs',LonSys='2',showbands=False,
                                coef=[0.,0.],subproj='',figxy=[8.0,4.0],
                                FiveMicron=True)
    
    @author: smhil
    """    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Data-Management-and-Access/processes')
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Data-Management-and-Access/Tests')

    import os
    import read_fits_V2 as RF2
    import numpy as np
    from astropy.io import fits
    sys.path.append('./Maps')
    #import read_fits_map_L2_L3 as RFM
    import make_patch as MP
    import map_and_context as mac
    import map_and_scatter as mas
    import map_and_scatter_SCubed as masc
    import map_cloudsurface as msurf
    import read_fits_V2 as RF2
    import L4_Jup_Map_Plot_V2 as L4M    
    ctbls=["terrain_r","Blues"]
    if  dataversion==1 or dataversion==2:
        fNH3low=60
        fNH3high=160
        PCldlow=1600
        PCldhigh=2200
    if dataversion=='H':
        fNH3low=0
        fNH3high=300
        PCldlow=1000
        PCldhigh=3000
        
    micronlow=0.5
    micronhigh=3.5
    print("############### len(obskey)= ",len(obskey))
    if len(obskey)==11:
        if (not FiveMicron) or FiveMicron=="png":
            if dataversion==2:
                PCldhdr,PClddata,fNH3hdr,fNH3data,RGB,RGB_CM,RGBtime= \
                                RF2.read_fits_map_L3_V2(obskey=obskey,
                                                        target=target,Level="L3",
                                                        LonSys=LonSys,
                                                        LimbCorrection=LimbCorrection,
                                                        dataversion=dataversion)
            elif dataversion==1:
                PCldhdr,PClddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime= \
                                RF2.read_fits_map_L3_V1(obskey=obskey,LonSys=LonSys,
                                                        Level="L3",
                                                        target=target,
                                                        LimbCorrection=LimbCorrection,
                                                        dataversion=dataversion)
            elif dataversion=="H":
                PCldhdr,PClddata,fNH3hdr,fNH3data,RGB,RGB_CM,RGBtime= \
                                RF2.read_fits_map_L3_V2(obskey=obskey,LonSys=LonSys,
                                                        Level="L3",
                                                        target=target,
                                                        LimbCorrection=LimbCorrection,
                                                        dataversion=dataversion)
            fNH3CM=fNH3hdr['CM'+LonSys]
            PCldCM=PCldhdr['CM'+LonSys]
        
    elif len(obskey)==17:
        print("@@@@@@@@@@@@@@@@")
        collection=obskey
        ###########################################################################
        # READ L4 BLENDED MAPS (fNH3, PCld, and IGB)
        fNH3data,fNH3stdv,fNH3frac,blendweightfNH3time,fNH3hdr,\
            PClddata,PCldstdv,PCldfrac,blendweightPCldtime,PCldhdr,\
            RGB,RGBstdv,RGBfrac,RGBtimearray,RGBhdr, \
            blendweightIRTF,IRTFstdv,IRTFfrac,IRTFtime,blendweightCH4889=\
            RF2.read_fits_map_L4(LonSys,collection=collection,
                                  collectionIRTF=False,
                                  collection889CH4=False)
        RGBtime=RGBhdr["DATE-OBS"]
        
    #elif FiveMicron=="fits":
    #    PCldhdr,PClddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime,micronhdr,microndatar= \
    #                    RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
    #                                            Level="L3",
    #                                            target=target,FiveMicron=FiveMicron)
                    
    pathmapplots='C:/Astronomy/Projects/SAS 2021 Ammonia/Data/L3 Plots/'+subproj+'/'
    if not os.path.exists(pathmapplots):
        os.makedirs(pathmapplots)
    ###########################################################################
    # Special for limb correction
    ###########################################################################             
    amfdata=1.0#(1.0/sza+1.0/eza)/2.0
    #figamf,axsamf=pl.subplots(figsize=(8.0,4.0), dpi=150, facecolor="white")
    #axsamf.imshow(amfdata,vmin=-5.,vmax=5.)
    if obskey=="20220730UTa":
        pathFITS='C:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
        amf=fits.open(pathFITS+"2022-07-30-amf_CM2_L360_MAP-BARE.fit")
        amf.info()
        amfhdr=amf[0].header
        amfdata=5.*amf[0].data/65535.
        amf.close()
    elif obskey=="20220919UTa":                
        pathFITS='C:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
        amf=fits.open(pathFITS+"2022-09-19-amf_CM2_L360_MAP-BARE.fit")
        amf.info()
        amfhdr=amf[0].header
        amfdata=5.*amf[0].data/65535.
        amf.close()

    #pl.imshow(amfdata)
    ###########################################################################
    # Set up figure and axes for plots
    ###########################################################################             
    if CMpref=='subobs':
        fNH3PlotCM=fNH3CM
        PCldPlotCM=PCldCM
        #if FiveMicron:
        #    micronPlotCM=micronhdr["CM"+LonSys]
    else:
        fNH3PlotCM=CMpref
        PCldPlotCM=CMpref
    NH3LonLims=[fNH3PlotCM-LonRng,fNH3PlotCM+LonRng]
    print("###################################################################")
    print("#######fNH3PlotCM=",fNH3PlotCM)
    print("fNH3PlotCM+LonRng,fNH3PlotCM-LonRng=",fNH3PlotCM+LonRng,fNH3PlotCM-LonRng)
    print("#######NH3LonLims=",NH3LonLims)
    print("#######360-NH3LonLims=",360-np.array(NH3LonLims))
    print("###################################################################")
    #amfpatch=MP.make_patch(amfdata,CoLatLims,NH3LonLims,fNH3PlotCM,LonRng,pad=True)

    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################
    #cbttl="Mean="+str(np.mean(fNH3_patch_mb))[:3]+" $\pm$ "+str(np.std(fNH3_patch_mb))[:2]
    fNH3_patch_mb,TestfNH3,tx_fNH3,fnNH3,RGB4Display=mac.map_and_context(fNH3data,
                                                       fNH3hdr["DATE-OBS"],fNH3hdr["BUNIT"],fNH3hdr["FILENAME"],
                                                       RGB,RGBtime,
                                                       LonSys,CoLatLims,NH3LonLims,
                                                       LonRng,fNH3PlotCM,
                                                       amfdata,coef[0],fNH3low,fNH3high,
                                                       showbands,FiveMicron,figxy,
                                                       ctbls[0],pathmapplots,Level='L3',cont=("contours" in plotoptions),
                                                       suptitle="Ammonia Mole Fraction",
                                                       cbar_title="Ammonia Mole Fraction (ppm)",
                                                       ROI=ROI,smoothcont=smoothcont,dataversion=dataversion)
    

    ###########################################################################
    ## Just RGB and Cloud Pressure
    ###########################################################################
    PCld_patch,TestPCld,tx_PCld,fnPCld,RGB4Display=mac.map_and_context(PClddata,
                                                        PCldhdr["DATE-OBS"],PCldhdr["BUNIT"],PCldhdr["FILENAME"],
                                                        RGB,RGBtime,
                                                        LonSys,CoLatLims,NH3LonLims,
                                                        LonRng,PCldPlotCM,
                                                        amfdata,coef[1],PCldlow,PCldhigh,
                                                        showbands,FiveMicron,figxy,
                                                        ctbls[1],pathmapplots,Level='L3',cont=("contours" in plotoptions),
                                                        suptitle="Cloud Top Pressure",
                                                        cbar_rev=True,
                                                        cbar_title="Cloud Top Pressure (mb)",
                                                        ROI=ROI,smoothcont=smoothcont,dataversion=dataversion)

    if "surface" in plotoptions:
        path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/"+subproj+"/"
        msurf.map_cloudsurface(PCld_patch,fNH3_patch_mb,RGB4Display,
                               PCldhdr,fNH3hdr,RGBtime,
                               LonSys,CoLatLims,[360-NH3LonLims[1],360-NH3LonLims[0]],
                               180,180,path,dataversion=dataversion)
        
    if "wave" in plotoptions:
        L4M.RossbyWavePlot(obskey,NH3LonLims,fNH3_patch_mb,PCld_patch,[6,6],
                           pathmapplots,LonSys,dataversion=dataversion)


    ###########################################################################
    ## Just RGB and 5 micron
    ###########################################################################
    if FiveMicron:
        micron_patch,Testmicron,tx_micron,fn5um=mac.map_and_context(np.log10(microndatar),micronhdr,
                                                          RGB,RGBtime,
                                                          LonSys,CoLatLims,NH3LonLims,
                                                          LonRng,PCldPlotCM,
                                                          amfdata,0.0,micronlow,micronhigh,
                                                          showbands,FiveMicron,figxy,
                                                          "gist_heat",pathmapplots,Level='L3',
                                                          suptitle="5 micron Radiance (log10)",
                                                          cbar_title="Log10(5um radiance)")

   
    ###########################################################################
    ## Compute Band or ROI Scatter Plot (PCloud vs fNH3)
    ###########################################################################
    #mas.map_and_scatter(fNH3_patch_mb,PCld_patch,PClddata,fNH3hdr['DATE-OBS'],LonSys,

    if "scatter" in plotoptions:
        dateobs,roilabel,mean1,stdv1,mean2,stdv2= \
            masc.map_and_scatter_SCubed(fNH3_patch_mb,PCld_patch,PClddata,RGB4Display,fNH3hdr['DATE-OBS'],LonSys,
            CoLatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
            amfdata,coef[0],tx_fNH3,fNH3low,fNH3high,PCldlow,PCldhigh,
            figxy,ctbls[1],pathmapplots,"PCloud & fNH3 (contours)",
            "PCloud vs fNH3",Level='L3',cbar_rev=True,cbar_title="Cloud-top Pressure (mb)",
            axis_inv=True,ROI=ROI,cont=("contours" in plotoptions),smoothcont=smoothcont,dataversion=dataversion)
        ROIout={obskey:{'dateobs':dateobs,'roilabel':roilabel,'mean1':mean1,'stdv1':stdv1,
                'mean2':mean2,'stdv2':stdv2}}#,'meanamf':meanamf}}
        print("############### ROIout= ",ROIout)
        
    if 'resid' in plotoptions:
        LonLimsEast=[360-NH3LonLims[1],360-NH3LonLims[0]]

        figrc,axsrc=pl.subplots(2,1,figsize=(5,5), dpi=150, facecolor="white",
                                sharex=True,sharey=True)
        figrc.suptitle(obskey+" Normalized Residuals")

        for i in [0,1]:
            axsrc[i].grid(linewidth=0.2)
            axsrc[i].ylim=[-45.,45.]
            axsrc[i].xlim=[360-LonLimsEast[0],360-LonLimsEast[1]]
            axsrc[i].set_xticks(np.linspace(450,0,31), minor=False)
            xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
            axsrc[i].set_xticklabels(xticklabels.astype(int))
            axsrc[i].set_yticks(np.linspace(-45,45,7), minor=False)
            axsrc[i].tick_params(axis='both', which='major', labelsize=9)
            axsrc[i].set_ylabel("PG Lat. (deg)",fontsize=10)
            #axsrc[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
        axsrc[0].set_title("Normalized Residual",fontsize=10)
        axsrc[1].set_title("Box Correlation (5x5 deg)",fontsize=10)


        norm_A,norm_B,resid_AB=residual_2d(fNH3_patch_mb,PCld_patch)
        if dataversion=='H':
            vn,vx=-0.6,0.6
        else:
            vn,vx=-1.0,1.0
        pp.plot_patch(resid_AB,CoLatLims,LonLimsEast,PCldPlotCM,LonRng,'BrBG',axsrc[0],
                   cbarplot=True,cbar_title="Test",cbar_reverse=False,vn=vn,vx=vx,n=6)
        
        if 'correl' in plotoptions:
            #figc,axsc=pl.subplots(1,figsize=(5,5), dpi=150, facecolor="white")
            #figc.suptitle(obskey+" Correlation")
    
            cc2d=rolling_corr2d(norm_A,norm_B,5,dataversion=dataversion)
            print("########### cc2d.shape= ",cc2d.shape,np.max(cc2d))
            pp.plot_patch(cc2d,CoLatLims,LonLimsEast,PCldPlotCM,LonRng,'seismic_r',axsrc[1],
                       cbarplot=True,cbar_title="Test",cbar_reverse=False,vn=-1.0,vx=1.0,n=6,alpha=0.7)
            
        if ROI:
            for R in ROI:
                for i in [0,1]:
                    axsrc[i].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                                  ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                                  ROI[R][2]+ROI[R][3]]),
                                  90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                                  ROI[R][1],ROI[R][0]]))
    ###########################################################################
    ## Compute Scatter Plot (PCloud vs 5um radiance)
    ###########################################################################
    if FiveMicron:
        mas.map_and_scatter(PCld_patch,micron_patch,np.log10(microndatar),fNH3hdr,LonSys,
                        CoLatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
                        coef[1],tx_PCld,PCldlow,PCldhigh,micronlow,micronhigh,
                        figxy,"gist_heat",pathmapplots,"5um Radiance & PCloud (contours)",
                        "PCloud vs 5um Radiance",Level='L3',FiveMicron=True,cbar_rev=False,swap_xy=True,
                        axis_inv=True,cbar_title="Log10(5um radiance)")

    ###########################################################################
    ## Compute Scatter Plot (fNH3 vs 5um radiance)
    ###########################################################################
    if FiveMicron:
        mas.map_and_scatter(fNH3_patch_mb,micron_patch,np.log10(microndatar),fNH3hdr,LonSys,
                        CoLatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
                        0.0,tx_fNH3,fNH3low,fNH3high,micronlow,micronhigh,
                        figxy,"gist_heat",pathmapplots,"fNH3 vs 5um Radiance",
                        "5um Radiance & fNH3 (contours)",Level='L3',FiveMicron=True,
                        cbar_rev=False,swap_xy=False,
                        axis_inv=True,cbar_title="Log10(5um radiance)")
   
    #return(fig1,axs1,fig2,axs2,fig3,axs3)
    if "scatter" in plotoptions:
        return(ROIout)
    #return(dateobs,roilabel,mean1,stdv1,mean2,stdv2,meanamf)


def load_png(file_path):
    """
    Purpose: Properly load a 48-bit PNG file
    Read from KITTI .png file
    Args:
        file_path string: file path(absolute)
    Returns:
        data (numpy.array): data of image in (Height, Width, 3) layout
    
    FROM: https://www.programcreek.com/python/example/98900/png.Reader
    """
    import png
    import numpy as np

    flow_object = png.Reader(filename=file_path)
    flow_direct = flow_object.asDirect()
    flow_data = list(flow_direct[2])
    (w, h) = flow_direct[3]['size']

    flow = np.zeros((h, w, 3), dtype=np.float64)
    for i in range(len(flow_data)):
        flow[i, :, 0] = flow_data[i][0::3]
        flow[i, :, 1] = flow_data[i][1::3]
        flow[i, :, 2] = flow_data[i][2::3]

    return flow.astype(np.uint16) 

