def map_and_scatter_SCubed(patchx,patchy,mapydata,RGBpatch,dateobs,LonSys,
                    LatLims,LonLimsWest,LonRng,PlotCM,fnout,
                    amfdata,coef,txin,xlow,xhigh,ylow,yhigh,figxy,
                    ct,pathout,Ltitle,Rtitle,Level='L3',cont=False,FiveMicron=False,
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
    ct : TYPE
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

    print("############################ dataversion= ",dataversion)
    LonLimsEast=[360-LonLimsWest[1],360-LonLimsWest[0]]

    ###########################################################################
    ## Compute Scatter Plot (PCloud vs fNH3)
    ###########################################################################
    #fig3,axs3=pl.subplots(1,2,figsize=(figxy[0],figxy[1]), dpi=150, facecolor="white")
    #fig3.suptitle(suptitle,x=0.5,ha='center',color='k')
    #fig3.suptitle(dateobs.replace("_"," ")+", CM"+LonSys+"="
    #              +str(int(PlotCM)),x=0.5,ha='center',color='k')

    #fig3 = pl.figure(layout="constrained",figsize=(8,6))
    #fig3 = pl.figure(figsize=(10,6),constrained_layout=False)
    #subfigs = fig3.subfigures(1, 2, wspace=0.05)
    #axs3 = subfigs[0].subplots(3, 1,sharex=True)
    
    #gs = fig3.add_gridspec(1, 2, left=0.02, right=0.98, wspace=0.02)

    #subfig_left  = fig3.add_subfigure(gs[0],left=0.02)
    #subfig_right = fig3.add_subfigure(gs[1])

    #axs3 = subfig_left.subplots(3, 1, sharex=True)
    
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
    #subfigs[0].set_facecolor('lightblue')
    #subfigs[0].suptitle('subfigs[0]\nLeft side')
    #subfigs[0].supxlabel('xlabel for subfigs[0]')
    
    #axs1 = subfigs[1].subplots(1)
    #axs1 = subfig_right.subplots(1)

    axs1.set_title("Ammonia Mole Fraction versus Cloud Pressure",fontsize=12)
    axs1.set_box_aspect(1)
    #subfigs[1].suptitle('subfigs[1]')
    #subfigs[1].supylabel('ylabel for subfigs[1]')

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
    axs3[0].set_title("Cloud Pressure (mb)",fontsize=10,y=1.0)

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
    axs3[1].set_title("Context Image (673/502/395nm)",fontsize=10,y=0.98)

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
    axs3[2].set_title("Ammonia Mole Fraction (ppm)",fontsize=10,y=0.95)

    #axs3[0].set_adjustable('box') 
    #axs3[1].set_adjustable('box') 

    #Pcloud_patch,vn,vx,tx=PP.plot_patch(PClouddata,LatLims,NH3LonLimsEast,
    #                                 PCldPlotCM,LonRng,"jet",
    #                                 axs2[0],'%3.2f',cont=False,
    #                                 cbar_reverse=True,vn=400,vx=900,n=6)
    #Testy_patch=MP.make_patch(mapydata,LatLims,LonLimsEast,PlotCM,LonRng)

    cbttl="Mean="+str(np.mean(patchy))[:4]+" $\pm$ "+str(np.std(patchy))[:3]
    tp,vn,vx,tx,cbary=PP.plot_patch(patchy,LatLims,LonLimsEast,
                                     PlotCM,LonRng,ct,
                                     axs3[0],'%3.2f',
                                     cbar_reverse=cbar_rev,vn=ylow,vx=yhigh,n=6,
                                     cbar_title=cbttl)
    #Testy_patch=MP.make_patch(mapydata,LatLims,LonLimsEast,PlotCM,LonRng)
    tp,vn,vx,tx,cbarRGB=PP.plot_patch(RGBpatch,LatLims,LonLimsEast,
                                     PlotCM,LonRng,"terrain_r",
                                     axs3[1],'%3.2f',
                                     cbar_reverse=False,vn=xlow,vx=xhigh,n=6,
                                     cbar_title=cbar_title,cbarvis=False)

    #Testy_patch=MP.make_patch(mapydata,LatLims,LonLimsEast,PlotCM,LonRng)
    cbttl="Mean="+str(np.mean(patchx))[:3]+" $\pm$ "+str(np.std(patchx))[:2]
    tp,vn,vx,tx,cbarx=PP.plot_patch(patchx,LatLims,LonLimsEast,
                                     PlotCM,LonRng,"terrain_r",
                                     axs3[2],'%3.2f',
                                     cbar_reverse=False,vn=xlow,vx=xhigh,n=6,
                                     cbar_title=cbttl)

    #for ax in axs3:
    #    ax.set_anchor('W')    
    #cbary.set_anchor('W')
    #cbarRGB.set_anchor('W')
    #cbarx.set_anchor('W')

    if cont:
        patchxsmth = gaussian_filter(patchx, sigma=smoothcont)
        temp=PC.plot_contours_on_patch(axs3[0],patchxsmth,LatLims,LonLimsEast,
                                       txin,frmt='%3.0f',clr='k')

    if coef==0.0:
        correction='_C0'
    else:
        correction='_C1'
    
    if Level=='L2':
        fnout=fnout.replace('TNH3',Rtitle)
    elif Level=='L3':
        fnout=fnout.replace('fNH3',Rtitle)
    print(patchx.shape,patchy.shape)
    
    if swap_xy and not ROI:
        roilabel,mean1,stdv1,mean2,stdv2,BZ=pms.plot_map_scatter(patchx,patchy,PlotCM,
                 LatLims,axs1,xlow,xhigh,ylow,yhigh,FiveMicron,axis_inv=axis_inv,dataversion=dataversion)
        print("Case 1")
    if not swap_xy and not ROI:     
        roilabel,mean1,stdv1,mean2,stdv2,BZ=pms.plot_map_scatter(patchy,patchx,PlotCM,
                 LatLims,axs1,ylow,yhigh,xlow,xhigh,FiveMicron,axis_inv=axis_inv,dataversion=dataversion)
        print("Case 2")
    
    #print("##################### mean1,stdv1,mean2,stdv2= ",mean1,stdv1,mean2,stdv2)
    #print("##################### np.mean(patchx),np.mean(patchy)= ",np.mean(patchx),np.mean(patchy))
    
    if swap_xy and ROI:
        print("Calling ROI")
        roilabel,mean1,stdv1,mean2,stdv2=prs.plot_roi_scatter(patchx,patchy,PlotCM,
                 LatLims,LonLimsEast,axs1,xlow,xhigh,ylow,yhigh,FiveMicron,
                 axis_inv=axis_inv,ROI=ROI,amfpatch=amfdata,dataversion=dataversion)
    if not swap_xy and ROI:    
        print("Calling ROI")
        roilabel,mean1,stdv1,mean2,stdv2=prs.plot_roi_scatter(patchy,patchx,PlotCM,
                 LatLims,LonLimsEast,axs1,ylow,yhigh,xlow,xhigh,FiveMicron,
                 axis_inv=axis_inv,ROI=ROI,amfpatch=amfdata,dataversion=dataversion)
        
        
    axs3[1].tick_params(axis='both', which='major', labelsize=9)

    ROIcolors={"Hot Spot":'r',
         "Gyre":'g',
         "Cloud Plume":'b',
         "NEB Reference":'k'}           
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
        #patch1=patch1*1000.
    
        #figcor,axscor=pl.subplots(1,1,figsize=(6.0,4.), dpi=150, facecolor="white",
        #                    sharey=True,sharex=True)          
    
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

    #box = axs3[1].get_position()
    #axs3[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 0.5, box.height * 1.015])    
    #fig3.subplots_adjust(left=0.02, right=0.95, top=0.92, bottom=0.10)
    #subfig_left.subplots_adjust(left=0.01, right=0.95, top=0.92, bottom=0.10)
    #for ax in axs3:
    #    ax.tick_params(axis='y', pad=1)
    for ax in axs3:
        print(ax.get_position())
    fig3.savefig(pathout+fnout[:-4]+' scatter.png',dpi=300)
    
    if not ROI:
        meanamf=False

    return(dateobs,roilabel,mean1,stdv1,mean2,stdv2)#,meanamf)