def plot_roi_scatter(obskey,dateobs,patch1,patch2,Real_CM2,LatLims,LonLims,axscor,PCldlow,PCldhigh,
                 fNH3low,fNH3high,FiveMicron,axis_inv=False,ROI=False,amfpatch=False,
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
    from scipy import stats
    from matplotlib.patches import Ellipse
    
    ###########################################################################
    # LOOP OVER ROIS AND PLOT SCATTER IN APPROPRIATE COLOR
    ###########################################################################
    mean1=[]
    stdv1=[]
    mean2=[]
    stdv2=[]
    meanamf=[]
    roilabel=[]
    slope=[]
    intercept=[]
    r_value=[]
    p_value=[]
    std_err=[]
    cov_matrix=[]
    if dataversion=="H":
        scale=20
    else:
        scale=1
    ROIcolors={"Hot Spot":'r',
         "Gyre":'g',
         "Cloud Plume":'b',
         "Reference":'k'}           

    if ROI:
        mean1.append(np.mean(patch1))
        mean2.append(np.mean(patch2))
        #meanamf.append(np.mean(amfsubpatch))
        stdv1.append(np.std(patch1))
        stdv2.append(np.std(patch2))
        roilabel.append('All')
        slopei, intercepti, r_valuei, p_valuei, std_erri = stats.linregress(patch1.ravel(), patch2.ravel())
        slope.append(slopei)
        intercept.append(intercepti)
        r_value.append(r_valuei)
        p_value.append(p_valuei)
        std_err.append(std_erri)
        
        cov_matrixi = np.cov([patch1.ravel(), patch2.ravel()],rowvar=True)
        cov_matrix.append(cov_matrixi)
        # Calculate eigenvalues and eigenvectors for the ellipse geometry
        # eigh is optimized for symmetric matrices like covariance matrices
        eigenvaluesi, eigenvectorsi = np.linalg.eigh(cov_matrixi)
        # Sort them in descending order so the largest eigenvalue dictates the major axis
        order = eigenvaluesi.argsort()[::-1]
        eigenvaluesi = eigenvaluesi[order]
        eigenvectorsi = eigenvectorsi[:, order]
        # Calculate dimensions and rotation angle
        # Dimensions represent 1 standard deviation along the principal axes
        widthi = 2 * np.sqrt(eigenvaluesi[1])
        heighti = 2 * np.sqrt(eigenvaluesi[0])
        # Angle of rotation in degrees (from the first eigenvector)
        anglei = np.degrees(np.arctan2(eigenvectorsi[1, 0], -eigenvectorsi[0, 0]))

        ellipse=Ellipse(
            xy=(np.mean(patch2), np.mean(patch1)),
            width=widthi * 1.96,
            height=heighti * 1.96,
            angle=anglei,
            edgecolor='grey',
            alpha=0.2,
            facecolor="none",
            linestyle="--",
            linewidth=1.5,
            label="95% confidence")
        axscor.add_patch(ellipse)        

        
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
            ###################################################################
            # NEW WORK FOR STATISTICAL ANALYSIS - SMH 6/9/2026
            ###################################################################
            # "ravel" to create 1D numpy arrays
            slopei, intercepti, r_valuei, p_valuei, std_erri = stats.linregress(subpatch1.ravel(), subpatch2.ravel())
            slope.append(slopei)
            intercept.append(intercepti)
            r_value.append(r_valuei)
            p_value.append(p_valuei)
            std_err.append(std_erri)
            
            cov_matrixi = np.cov([subpatch1.ravel(), subpatch2.ravel()],rowvar=True)
            cov_matrix.append(cov_matrixi)
            # Calculate eigenvalues and eigenvectors for the ellipse geometry
            # eigh is optimized for symmetric matrices like covariance matrices
            eigenvaluesi, eigenvectorsi = np.linalg.eigh(cov_matrixi)
            # Sort them in descending order so the largest eigenvalue dictates the major axis
            order = eigenvaluesi.argsort()[::-1]
            eigenvaluesi = eigenvaluesi[order]
            eigenvectorsi = eigenvectorsi[:, order]
            # Calculate dimensions and rotation angle
            # Dimensions represent 1 standard deviation along the principal axes
            widthi = 2 * np.sqrt(eigenvaluesi[1])
            heighti = 2 * np.sqrt(eigenvaluesi[0])
            # Angle of rotation in degrees (from the first eigenvector)
            anglei = np.degrees(np.arctan2(eigenvectorsi[1, 0], -eigenvectorsi[0, 0]))

            print("################")
            print(np.mean(subpatch1), np.mean(subpatch2), widthi, heighti, anglei)
            
            ellipse=Ellipse(
                xy=(np.mean(subpatch2), np.mean(subpatch1)),
                width=widthi * 1.96,
                height=heighti * 1.96,
                angle=anglei,
                edgecolor=clr,
                facecolor="none",
                linestyle="--",
                linewidth=1.5,
                label="95% confidence")
            axscor.add_patch(ellipse)        
     
    axscor.grid(linewidth=0.2)
    axscor.set_ylim(PCldlow,PCldhigh)
    axscor.set_xlim(fNH3low,fNH3high)
    #axscor.set_ylabel("Cloud-top Pressure (mb)",fontsize=10)
    axscor.set_xlabel(xaxistitle,fontsize=10)
    axscor.set_ylabel(yaxistitle,fontsize=10)
    if axis_inv:
        axscor.invert_yaxis()
    """
    if FiveMicron:
        axscor.set_xlabel("5um Radiance (Log10(arb. units)",fontsize=10)
    else:
        axscor.set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=10)
    """
                    
    axscor.legend(fontsize=8,ncols=2,labelcolor='mfc')
    
    ROIout={obskey:{'dateobs':dateobs,'roilabel':roilabel,'mean1':mean1,'stdv1':stdv1,
            'mean2':mean2,'stdv2':stdv2,'slope':slope,'intercept':intercept,
            'r_value':r_value,'p_value':p_value,'std_err':std_err,'cov_matrix':cov_matrix}}#,'meanamf':meanamf}}

    return ROIout
    #return(roilabel,mean1,stdv1,mean2,stdv2,slope,intercept,r_value,p_value,std_err,cov_matrix)#,meanamf)
  
