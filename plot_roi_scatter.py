def plot_roi_scatter(patch1,patch2,Real_CM2,LatLims,LonLims,axscor,PCldlow,PCldhigh,
                 fNH3low,fNH3high,FiveMicron,axis_inv=False,ROI=False,amfpatch=False,
                 dataversion=2):
    """
    PURPOSE:    Takes two map patches and makes a scatter plot
    CALLS:      n/a
    CALLED BY:  map_and_scatter

    Parameters
    ----------
    patch1 : TYPE
        DESCRIPTION.
    patch2 : TYPE
        DESCRIPTION.
    Real_CM2 : TYPE
        DESCRIPTION.
    LatLims : TYPE
        DESCRIPTION.
    LonLims : TYPE
        DESCRIPTION.
    axscor : TYPE
        DESCRIPTION.
    PCldlow : TYPE
        DESCRIPTION.
    PCldhigh : TYPE
        DESCRIPTION.
    fNH3low : TYPE
        DESCRIPTION.
    fNH3high : TYPE
        DESCRIPTION.
    FiveMicron : TYPE
        DESCRIPTION.
    axis_inv : TYPE, optional
        DESCRIPTION. The default is False.
    ROI : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    import pylab as pl
    import numpy as np
    import copy
    
    ###########################################################################
    # LOOP OVER ROIS AND PLOT SCATTER IN APPROPRIATE COLOR
    ###########################################################################
    mean1=[]
    stdv1=[]
    mean2=[]
    stdv2=[]
    meanamf=[]
    roilabel=[]
    if dataversion=="H":
        scale=20
    else:
        scale=1
    ROIcolors={"Hot Spot":'r',
         "Gyre":'g',
         "Cloud Plume":'b',
         "NEB Reference":'k'}           

    if ROI:
        for R in ROI:
            RLatLims=-LatLims[0]+np.array([ROI[R][0],ROI[R][1]])
            RCM=ROI[R][2]
            RLonRng=ROI[R][3]
            RLonLims=[360-int(RCM+RLonRng),360-int(RCM-RLonRng)]

            RLonLims=np.array(RLonLims)-LonLims[0]

            subpatch1=patch1[RLatLims[0]*scale:RLatLims[1]*scale,
                             RLonLims[0]*scale:RLonLims[1]*scale]
            subpatch2=patch2[RLatLims[0]*scale:RLatLims[1]*scale,
                             RLonLims[0]*scale:RLonLims[1]*scale]
            #amfsubpatch=amfpatch[RLatLims[0]*scale:RLatLims[1]*scale,
            #                 RLonLims[0]*scale:RLonLims[1]*scale]

            clr=ROIcolors[R]
            if dataversion=="H":
                axscor.scatter(subpatch2,subpatch1,marker=".",s=2,color=clr,linewidths=0,alpha=1.0,label=R)
            else:
                axscor.scatter(subpatch2,subpatch1,marker="o",s=3.0,color=clr,alpha=0.8,label=R)
            
            mean1.append(np.mean(subpatch1))
            mean2.append(np.mean(subpatch2))
            #meanamf.append(np.mean(amfsubpatch))
            stdv1.append(np.std(subpatch1))
            stdv2.append(np.std(subpatch2))
            roilabel.append(R)
     
    axscor.grid(linewidth=0.2)
    axscor.set_ylim(PCldlow,PCldhigh)
    axscor.set_xlim(fNH3low,fNH3high)
    axscor.set_ylabel("Cloud-top Pressure (mb)",fontsize=10)
    if axis_inv:
        axscor.invert_yaxis()
    if FiveMicron:
        axscor.set_xlabel("5um Radiance (Log10(arb. units)",fontsize=10)
    else:
        axscor.set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=10)
                    
    axscor.legend(fontsize=8,ncols=2,labelcolor='mfc')
    
    
    return(roilabel,mean1,stdv1,mean2,stdv2)#,meanamf)
  
