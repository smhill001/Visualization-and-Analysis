def plot_map_scatter(obskey,dateobs,patch1,patch2,Real_CM2,LatLims,axscor,PCldlow,PCldhigh,
                 fNH3low,fNH3high,FiveMicron,axis_inv=False,Bands=False,
                 dataversion=2,xaxistitle='',yaxistitle=''):
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
    Bands : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    import copy
    import plot_roi_scatter as prs
    
    #print("000000000: ",patch1.shape,patch2.shape)
    ###########################################################################
    # SET BELT AND ZONE BOUNDARIES
    ###########################################################################
    BZ={"SSTB":[-39.6,-36.2],
          "STZ":[-36.2,-32.4],
          "STB":[-32.4,-27.1],
          "STrZ":[-27.1,-19.7],
          "SEB":[-19.7,-7.2],
          "SEZ":[-7.2,0.0],
          "NEZ":[0.0,6.9],
          "NEB":[6.9,17.4],
          "NTrZ":[17.4,24.2],
          "NTB":[24.2,31.4],
          "NTZ":[31.4,35.4],
          "NNTB":[35.4,39.6]}

    BZind=copy.deepcopy(BZ)   
    BZkeys=BZ.keys()

    if dataversion=="H":
        scale=20
    else:
        scale=1

    ###########################################################################
    # LOOP OVER BELTS AND PLOT SCATTER IN APPROPRIATE COLOR
    ###########################################################################   
    ROI={}
    
    ROIout={obskey:{'dateobs':dateobs,'roilabel':[],'nsamples':[],'mean1':[],'stdv1':[],
            'mean2':[],'stdv2':[],'slope':[],'intercept':[],
            'r_value':[],'p_value':[],'std_err':[],'cov_matrix':[]}}

    counter=0
    for key in BZ.keys():
        #print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
        BZind[key][0]=int(90-BZ[key][0])-LatLims[0]
        BZind[key][1]=int(90-BZ[key][1])-LatLims[0]
        
        if BZind[key][0]<0:
            BZind[key][0]=0
        if BZind[key][1]<0:
            BZind[key][1]=0
        if BZind[key][0]>(LatLims[1]-LatLims[0]):
            BZind[key][0]=LatLims[1]-LatLims[0]
        if BZind[key][1]>(LatLims[1]-LatLims[0]):
            BZind[key][1]=LatLims[1]-LatLims[0]
        
        if BZind[key][0]==BZind[key][1]:
            print("do nothing")
        else:
            counter=counter+1
            ROI[key]=[BZind[key][1],BZind[key][0]]
            ROIcolor='C'+str(counter)
            print(patch2.shape,patch1.shape)
            subpatch1=patch1[BZind[key][1]*scale:BZind[key][0]*scale,:]
            subpatch2=patch2[BZind[key][1]*scale:BZind[key][0]*scale,:]
            if dataversion=='H':
                axscor.scatter(subpatch2,subpatch1,
                               marker=".",s=0.1,linewidths=0,
                               alpha=1.0,label=key)

            else:
                axscor.scatter(subpatch2,subpatch1,
                               marker="o",s=3.0,
                               alpha=0.8,label=key)
            R=[BZind[key][0],BZind[key][1],]
            ROIout=prs.statistics_helper(ROIout,obskey,R,subpatch1,subpatch2,
                                     ROIcolor,axscor,alpha=1.0)

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
        
    axscor.set_xlabel=xaxistitle
    axscor.set_ylabel=yaxistitle
                    
    axscor.legend(fontsize=7,ncols=4,labelcolor='mfc')
    
    #return(keylabel,mean1,stdv1,mean2,stdv2,BZ)
    return ROIout,BZ
  
