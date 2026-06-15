def map_and_scatter_SCubed(obskey,patchx,patchy,mapydata,RGBpatch,dateobs,LonSys,
                    LatLims,LonLimsWest,LonRng,PlotCM,fnout,
                    amfdata,coef,txinx,txiny,xlow,xhigh,ylow,yhigh,figxy,
                    ctbls,pathout,Ltitle,Rtitle,Level='L3',
                    maptitles=['','',''],
                    cont=False,FiveMicron=False,
                    cbar_rev=False,swap_xy=False,axis_inv=False,cbar_title="Test",
                    suptitle="Test",ROI=False,smoothcont=0,dataversion=2):
    """
    PURPOSE:    Makes a pair of plots, the left one is a patch map of one data
                set overlayed by another patch map data set. The right plot is
                a scatter plot of the two patch maps.
                
    CALLS:      plot_map_scatter
                plot_roi_scatter

    Parameters
    ----------
    patchx : TYPE
        DESCRIPTION.
    patchy : TYPE
        DESCRIPTION.
    mapydata : TYPE
        DESCRIPTION.
    dateobs : TYPE
        DESCRIPTION.
    LonSys : TYPE
        DESCRIPTION.
    LatLims : TYPE
        DESCRIPTION.
    LonLims : TYPE
        DESCRIPTION.
    LonRng : TYPE
        DESCRIPTION.
    PlotCM : TYPE
        DESCRIPTION.
    fnout : TYPE
        DESCRIPTION.
    coef : TYPE
        DESCRIPTION.
    txin : TYPE
        DESCRIPTION.
    xlow : TYPE
        DESCRIPTION.
    xhigh : TYPE
        DESCRIPTION.
    ylow : TYPE
        DESCRIPTION.
    yhigh : TYPE
        DESCRIPTION.
    figxy : TYPE
        DESCRIPTION.
    ctbls : TYPE
        DESCRIPTION.
    pathout : TYPE
        DESCRIPTION.
    Ltitle : TYPE
        DESCRIPTION.
    Rtitle : TYPE
        DESCRIPTION.
    Level : TYPE, optional
        DESCRIPTION. The default is 'L3'.
    FiveMicron : TYPE, optional
        DESCRIPTION. The default is False.
    cbar_rev : TYPE, optional
        DESCRIPTION. The default is False.
    swap_xy : TYPE, optional
        DESCRIPTION. The default is False.
    axis_inv : TYPE, optional
        DESCRIPTION. The default is False.
    cbar_title : TYPE, optional
        DESCRIPTION. The default is "Test".
    suptitle : TYPE, optional
        DESCRIPTION. The default is "Test".
    ROI : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import pylab as pl
    import numpy as np
    sys.path.append('./Maps')
    import plot_patch as PP
    import make_patch as MP
    import plot_contours_on_patch as PC
    import copy
    import plot_map_scatter as pms
    import plot_roi_scatter as prs
    from scipy.ndimage import gaussian_filter

    statistics={}
    print("############################ dataversion= ",dataversion)
    LonLimsEast=[360-LonLimsWest[1],360-LonLimsWest[0]]

    ###########################################################################
    ## Compute Scatter Plot (PCloud vs fNH3)
    ###########################################################################
  
    fig3 = pl.figure(figsize=(8,4.5),dpi=150,facecolor="white")


    gs = fig3.add_gridspec(
        3, 2,
        left=-0.05,
        right=0.98,
        top=0.94,
        bottom=0.10,
        wspace=0.0,
        hspace=0.4
    )

    # Left column (your 3 stacked plots)
    axs3 = [fig3.add_subplot(gs[i, 0]) for i in range(3)]
    axs3[0].sharex(axs3[2])
    axs3[2].sharex(axs3[2])
    axs3[0].sharey(axs3[2])
    axs3[2].sharey(axs3[2])
    
    # Right column (whatever goes there)
    axs1 = fig3.add_subplot(gs[:, 1])
    axs1.set_title("Ammonia Mole Fraction versus Cloud Pressure",fontsize=12)
    axs1.set_box_aspect(1)

    axs3[0].grid(linewidth=0.2)
    axs3[0].ylim=[-45.,45.]
    axs3[0].xlim=[360-LonLimsEast[0],360-LonLimsEast[1]]
    axs3[0].set_xticks(np.linspace(450,0,31), minor=False)
    xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
    axs3[0].set_xticklabels(xticklabels.astype(int))
    axs3[0].set_yticks(np.linspace(-45,45,7), minor=False)
    axs3[0].tick_params(axis='both', which='major', labelsize=9)
    axs3[0].set_ylabel("PG Lat. (deg)",fontsize=10)
    #axs3[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs3[0].set_title(maptitles[0],fontsize=10,y=1.0)

    axs3[1].grid(linewidth=0.2)
    axs3[1].ylim=[-45.,45.]
    axs3[1].xlim=[360-LonLimsEast[0],360-LonLimsEast[1]]
    axs3[1].set_xticks(np.linspace(450,0,31), minor=False)
    xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
    axs3[1].set_xticklabels(xticklabels.astype(int))
    axs3[1].set_yticks(np.linspace(-45,45,7), minor=False)
    axs3[1].tick_params(axis='both', which='major', labelsize=9)
    axs3[1].set_ylabel("PG Lat. (deg)",fontsize=10)
    #axs3[1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs3[1].set_title(maptitles[1],fontsize=10,y=0.98)

    axs3[2].grid(linewidth=0.2)
    axs3[2].ylim=[-45.,45.]
    axs3[2].xlim=[360-LonLimsEast[0],360-LonLimsEast[1]]
    axs3[2].set_xticks(np.linspace(450,0,31), minor=False)
    xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
    axs3[2].set_xticklabels(xticklabels.astype(int))
    axs3[2].set_yticks(np.linspace(-45,45,7), minor=False)
    axs3[2].tick_params(axis='both', which='major', labelsize=9)
    axs3[2].set_ylabel("PG Lat. (deg)",fontsize=10)
    axs3[2].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs3[2].set_title(maptitles[2],fontsize=10,y=0.95)


    cbttl="Mean="+str(np.mean(patchy))[:4]+" $\pm$ "+str(np.std(patchy))[:3]
    statistics |={'mean_y':np.mean(patchy),'mean_x':np.mean(patchx),
                  'stdv_y':np.std(patchy),'stdv_x':np.std(patchx)}
    #cbttl='test while mean not working for CI or AOI'
    tp,vn,vx,tx,cbary=PP.plot_patch(patchy,LatLims,LonLimsEast,
                                     PlotCM,LonRng,ctbls[1],
                                     axs3[0],'%3.2f',
                                     cbar_reverse=cbar_rev,vn=ylow,vx=yhigh,n=6,
                                     cbar_title=cbttl)
    #Testy_patch=MP.make_patch(mapydata,LatLims,LonLimsEast,PlotCM,LonRng)
    tp,vn,vx,tx,cbarRGB=PP.plot_patch(RGBpatch,LatLims,LonLimsEast,
                                     PlotCM,LonRng,ctbls[0],
                                     axs3[1],'%3.2f',
                                     cbar_reverse=False,vn=xlow,vx=xhigh,n=6,
                                     cbar_title=cbar_title,cbarvis=False)

    #Testy_patch=MP.make_patch(mapydata,LatLims,LonLimsEast,PlotCM,LonRng)
    cbttl="Mean="+str(np.mean(patchx))[:3]+" $\pm$ "+str(np.std(patchx))[:2]
    tp,vn,vx,tx,cbarx=PP.plot_patch(patchx,LatLims,LonLimsEast,
                                     PlotCM,LonRng,ctbls[0],
                                     axs3[2],'%3.2f',
                                     cbar_reverse=False,vn=xlow,vx=xhigh,n=6,
                                     cbar_title=cbttl)

    if cont:
        patchxsmth = gaussian_filter(patchx, sigma=smoothcont)
        temp=PC.plot_contours_on_patch(axs3[2],patchxsmth,LatLims,LonLimsEast,
                                       txinx,frmt='%3.0f',clr='k')
        patchysmth = gaussian_filter(patchy, sigma=smoothcont)
        temp=PC.plot_contours_on_patch(axs3[0],patchysmth,LatLims,LonLimsEast,
                                       txiny,frmt='%3.0f',clr='r')

    if coef==0.0:
        correction='_C0'
    else:
        correction='_C1'
    
    if dataversion==2:
        if Level=='L2':
            fnout=fnout.replace('TNH3',Rtitle)
        elif Level=='L3':
            fnout=fnout.replace('fNH3',Rtitle)
        axs1.scatter(patchx,patchy,marker="o",s=3.0,color='grey',alpha=0.2,label='All')
    elif dataversion=='H':
        fnout=fnout.replace('.png',Rtitle+'.png')
        axs1.scatter(patchx,patchy,marker=".",s=2,color='grey',linewidths=0,alpha=0.2,label='All')
    
    print(patchx.shape,patchy.shape)
    
    
    if swap_xy and not ROI:
        roilabel,mean1,stdv1,mean2,stdv2,BZ=pms.plot_map_scatter(patchx,patchy,PlotCM,
                 LatLims,axs1,xlow,xhigh,ylow,yhigh,FiveMicron,axis_inv=axis_inv,
                 dataversion=dataversion,xaxistitle=maptitles[2],yaxistitle=maptitles[0])
        print("Case 1")
    if not swap_xy and not ROI:     
        roilabel,mean1,stdv1,mean2,stdv2,BZ=pms.plot_map_scatter(patchy,patchx,PlotCM,
                 LatLims,axs1,ylow,yhigh,xlow,xhigh,FiveMicron,axis_inv=axis_inv,
                 dataversion=dataversion,xaxistitle=maptitles[2],yaxistitle=maptitles[0])
        print("Case 2")
    
    #print("##################### mean1,stdv1,mean2,stdv2= ",mean1,stdv1,mean2,stdv2)
    #print("##################### np.mean(patchx),np.mean(patchy)= ",np.mean(patchx),np.mean(patchy))
    
    if swap_xy and ROI:
        print("Calling ROI",maptitles[2],maptitles[0])
        ROIout=prs.plot_roi_scatter(obskey,dateobs,patchx,patchy,PlotCM,
                 LatLims,LonLimsEast,axs1,xlow,xhigh,ylow,yhigh,FiveMicron,
                 axis_inv=axis_inv,ROI=ROI,amfpatch=amfdata,
                 dataversion=dataversion,xaxistitle=maptitles[2],yaxistitle=maptitles[0])
    if not swap_xy and ROI:    
        print("Calling ROI",maptitles[2],maptitles[0])
        ROIout=prs.plot_roi_scatter(obskey,dateobs,patchy,patchx,PlotCM,
                 LatLims,LonLimsEast,axs1,ylow,yhigh,xlow,xhigh,FiveMicron,
                 axis_inv=axis_inv,ROI=ROI,amfpatch=amfdata,
                 dataversion=dataversion,xaxistitle=maptitles[2],yaxistitle=maptitles[0])
        
        
    axs3[1].tick_params(axis='both', which='major', labelsize=9)

    ROIcolors={"Hot Spot":'r',
         "Gyre":'g',
         "Cloud Plume":'b',
         "Reference":'k'}           
    if ROI:
        for R in ROI:
            
            for i in [0,1,2]:
                axs3[i].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                              ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                              ROI[R][2]+ROI[R][3]]),
                              90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                              ROI[R][1],ROI[R][0]]),color=ROIcolors[R])
    else:
        BZind=copy.deepcopy(BZ)   
        BZkeys=BZ.keys()
        BZind=copy.deepcopy(BZ)   
        BZkeys=BZ.keys()
    
        clrind=0
        for key in BZ.keys():
            #print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
            #######################################################################
            # Compute the axis fraction and plot vertical bars using axvspan
            BZind[key][0]=1.-((90-BZ[key][0])-LatLims[0])/(LatLims[1]-LatLims[0])
            BZind[key][1]=1.-((90-BZ[key][1])-LatLims[0])/(LatLims[1]-LatLims[0])
            #print(key,BZind[key])
                
            #print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
            clr='C'+str(clrind)
            #print(clr)
            if BZind[key][0]>1.0 or BZind[key][1]<0.0:
                print("do nothing")
            else:
                axs3[0].axvspan(360-LonLimsEast[0],360-LonLimsEast[0]-1,
                                ymin=BZind[key][1],ymax=BZind[key][0],alpha=1.0,
                                color=clr)
                axs3[0].axvspan(360-LonLimsEast[1],360-LonLimsEast[1]+1,
                                ymin=BZind[key][1],ymax=BZind[key][0],alpha=1.0,
                                color=clr)
                
                clrind=clrind+1

    #fig3.savefig(pathout+fnout[:-4]+' scatter.png',dpi=300)
    
    return ROIout,fig3,axs1,axs3,fnout