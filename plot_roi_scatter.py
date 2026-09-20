import numpy as np
from scipy import stats
from matplotlib.patches import Ellipse
import matplotlib.lines as mlines



def pooled_covariance_ROIs(ROIout, obskey):

    covs = ROIout[obskey]['cov_matrix'][1:]   # exclude parent
    ns   = ROIout[obskey]['nsamples'][1:]

    num = np.zeros((2,2))
    den = 0

    for n, S in zip(ns, covs):
        num += (n - 1) * np.asarray(S)
        den += (n - 1)

    return num / den

def roi_pairwise_mahalanobis(ROIout, obskey):

    labels = ROIout[obskey]['roilabel'][1:]   # exclude parent

    means = np.column_stack([
        ROIout[obskey]['mean1'][1:],
        ROIout[obskey]['mean2'][1:]
    ])

    Sp = pooled_covariance_ROIs(ROIout, obskey)
    Sinv = np.linalg.inv(Sp)

    n = len(labels)
    D = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):

            d = means[i] - means[j]
            D2 = d.T @ Sinv @ d

            D[i, j] = np.sqrt(D2)
            D[j, i] = D[i, j]
            
    results={'labels':labels,'D':D}

    return results

def mahalanobis_to_parent(ROIout, obskey):

    means = np.column_stack([
        ROIout[obskey]['mean1'],
        ROIout[obskey]['mean2']
    ])

    mu0 = means[0]

    S0 = np.array(
        ROIout[obskey]['cov_matrix'][0]
    )

    Sinv = np.linalg.inv(S0)

    results = {}

    for i in range(1, len(means)):

        d = means[i] - mu0

        D2 = float(d.T @ Sinv @ d)

        results[
            ROIout[obskey]['roilabel'][i]
        ] = {
            'D2': D2,
            'D': np.sqrt(D2)
        }

    return results

def pairwise_mahalanobis(ROIout, obskey):

    means = np.column_stack([
        ROIout[obskey]['mean1'],
        ROIout[obskey]['mean2']
    ])

    labels = ROIout[obskey]['roilabel']

    S0 = np.array(
        ROIout[obskey]['cov_matrix'][0]
    )

    Sinv = np.linalg.inv(S0)

    n = len(labels)

    D = np.zeros((n, n))
    #!!!!! Should make the output an array of dicts like mahalanobis_to_parent
    for i in range(n):
        for j in range(i + 1, n):

            d = means[i] - means[j]

            D2 = d.T @ Sinv @ d

            D[i, j] = np.sqrt(D2)
            D[j, i] = D[i, j]
            
    results={'labels':labels,'D':D}

    return results

def plot_Mahal_ellipse(cov_matrixi,mean_patch1,mean_patch2,axscor,clr,alpha=0.8):
    # Calculate eigenvalues and eigenvectors for the ellipse geometry
    # eigh is optimized for symmetric matrices like covariance matrices
    eigenvaluesi, eigenvectorsi = np.linalg.eigh(cov_matrixi)
    # Sort them in descending order so the largest eigenvalue dictates the major axis
    order = eigenvaluesi.argsort()[::-1]
    eigenvaluesi = eigenvaluesi[order]
    eigenvectorsi = eigenvectorsi[:, order]
    # Calculate dimensions and rotation angle
    # Dimensions represent 1 standard deviation along the principal axes
    """
    widthi = 2 * np.sqrt(eigenvaluesi[1])
    heighti = 2 * np.sqrt(eigenvaluesi[0])
    # Angle of rotation in degrees (from the first eigenvector)
    anglei = np.degrees(np.arctan2(eigenvectorsi[1, 0], -eigenvectorsi[0, 0]))
    """
    # 2. Correct 2D 95% scale factor (~2.4477)
    from scipy.stats import chi2
    scale_factor = np.sqrt(chi2.ppf(0.95, df=2))
    
    # 3. Assign Major Axis to Width, Minor Axis to Height
    widthi  = 2 * np.sqrt(eigenvaluesi[0]) * scale_factor  # Major axis
    heighti = 2 * np.sqrt(eigenvaluesi[1]) * scale_factor  # Minor axis
    
    # 4. Correct arctan2(y, x) vector alignment
    # Row 0 = patch1 (Y-axis), Row 1 = patch2 (X-axis)
    y_component = eigenvectorsi[0, 0]
    x_component = eigenvectorsi[1, 0]
    anglei = np.degrees(np.arctan2(y_component, x_component))
    ellipse=Ellipse(
        xy=(mean_patch2, mean_patch1),
        #width=widthi * 1.96,
        #height=heighti * 1.96,
        width=widthi,
        height=heighti,
        angle=anglei,
        edgecolor=clr,
        alpha=alpha,
        facecolor="none",
        linestyle="solid",
        linewidth=1.0)#,
        #label="95% Confidence")
    axscor.add_patch(ellipse)
    ## 2. Create a clean Line2D proxy matching your ellipse styling
    #axscor.ellipse_proxy = mlines.Line2D(
    #    [], [], color=clr, linestyle="--", linewidth=1.5, label="95% Confidence")   
    
def statistics_helper(ROIout,obskey,R,patch1,patch2,clr,axscor,alpha=1.0,plot_ellipse=True):
    mean_patch1,mean_patch2,stdv_patch1,stdv_patch2 = np.mean(patch1), np.mean(patch2),np.std(patch1),np.std(patch2)
    ROIout[obskey]['mean1'].append(mean_patch1)
    ROIout[obskey]['mean2'].append(mean_patch2)
    ROIout[obskey]['stdv1'].append(stdv_patch1)
    ROIout[obskey]['stdv2'].append(stdv_patch2)
    ROIout[obskey]['roilabel'].append(R)
    ROIout[obskey]['nsamples'].append(patch1.size)
    
    #slopeiold, interceptiold, r_valueiold, p_valueiold, std_erriold = stats.linregress(patch1.ravel(), patch2.ravel())
    #print("###############")
    #print("slopeiold, interceptiold, r_valueiold, p_valueiold, std_erriold=",slopeiold, interceptiold, r_valueiold, p_valueiold, std_erriold)
    cov_matrixi = np.cov([patch1.ravel(), patch2.ravel()],rowvar=True)
    
    # 2. New code from Gemini for slope, intercept etc. from covariance matrix
    var_x, var_y = cov_matrixi[0, 0], cov_matrixi[1, 1]
    cov_xy = cov_matrixi[0, 1]
    # 3. Major Axis Slope & Intercept (Aligns with ellipse long axis)
    slopei = ((var_y - var_x) + np.sqrt((var_y - var_x)**2 + 4 * cov_xy**2)) / (2 * cov_xy)
    intercepti = mean_patch2 - slopei * mean_patch1
    
    # 4. Correlation r and p-value
    r_valuei = cov_xy / np.sqrt(var_x * var_y)
    t_stat = r_valuei * np.sqrt(patch1.size - 2) / np.sqrt(1 - r_valuei**2)
    p_valuei = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=patch1.size - 2))
    
    # 5. Standard error of slope via eigenvalues
    evals, evecs = np.linalg.eigh(cov_matrixi)
    l1, l2 = evals[1], evals[0]  # l1 is largest eigenvalue
    std_erri = (1 + slopei**2) * np.sqrt((l1 * l2) / (patch1.size * (l1 - l2)**2))
    print("slopei, intercepti, r_valuei, p_valuei, std_erri=",slopei, intercepti, r_valuei, p_valuei, std_erri)
    print("###############")
    # 3. Sort descending to get major axis
    #v_major = evecs[:, np.argmax(evals)]  # [fNH3_comp, PCld_comp]
    
    # 4. Physical slope = delta(Y) / delta(X) = delta(fNH3) / delta(PCld)
    #slope_physical = v_major[0] / v_major[1]  # Returns ppm / mb
    #print(slope_physical)
    #slope_sma = np.sign(cov_matrixi[0, 1]) * np.sqrt(cov_matrixi[0, 0] / cov_matrixi[1, 1])
    #print(slope_sma)
    slope_ols_ppm_mb = cov_xy / var_y
    intercept_ols_ppm = mean_patch1 - (slope_ols_ppm_mb * mean_patch2)
    print(slope_ols_ppm_mb,intercept_ols_ppm)
    df = patch1.size - 2
    std_err_ppm_mb = abs(slope_ols_ppm_mb) * np.sqrt((1 - r_valuei**2) / (df * r_valuei**2))
    ROIout[obskey]['slope'].append(slope_ols_ppm_mb)
    ROIout[obskey]['intercept'].append(intercept_ols_ppm)
    ROIout[obskey]['r_value'].append(r_valuei)
    ROIout[obskey]['p_value'].append(p_valuei)
    ROIout[obskey]['std_err'].append(std_err_ppm_mb)
    ROIout[obskey]['cov_matrix'].append(cov_matrixi)
    
    #print("cov_matrixi",cov_matrixi)
    if plot_ellipse:
        plot_Mahal_ellipse(cov_matrixi,mean_patch1,mean_patch2,axscor,clr,alpha=0.8)
    #return
    return ROIout


    
def plot_roi_scatter(obskey,dateobs,ROI_ID,patch1,patch2,Real_CM2,LatLims,LonLims,axscor,PCldlow,PCldhigh,
                 fNH3low,fNH3high,FiveMicron,axis_inv=False,ROI=False,amfpatch=False,
                 dataversion=2,xaxistitle='',yaxistitle='',fNH3factor=1.0,
                 plot_ellipse=True,GMM_clusters=0):
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
    import numpy as np
    print("plot_roi_scatter plot_roi_scatter ")
    print("DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG ")
    print("fNH3factor= ",fNH3factor)
    print("END END END END END END END END")

    ###########################################################################
    # Initial setup and statistics on whole population
    ###########################################################################

    if dataversion=="H":
        scale=20
    else:
        scale=1
    ROIcolors={"Hot Spot":'r',
         "Gyre":'g',
         "Cloud Plume":'b',
         "Reference":'k'}           

    ROIout={obskey:{'dateobs':dateobs,'roilabel':[],'nsamples':[],'mean1':[],'stdv1':[],
            'mean2':[],'stdv2':[],'slope':[],'intercept':[],
            'r_value':[],'p_value':[],'std_err':[],'cov_matrix':[]}}

    ROIout=statistics_helper(ROIout,obskey,'All',patch1,patch2*fNH3factor,
                             'grey',axscor,alpha=0.2,plot_ellipse=plot_ellipse)
    
    ftemp=np.linspace(0, 450,num=50)
    ptemp=np.array(ROIout[obskey]['intercept'])+np.array(ROIout[obskey]['slope'])*np.array(ftemp)
    print("ROIout[obskey]['intercept'],ROIout[obskey]['slope']=",ROIout[obskey]['intercept'],ROIout[obskey]['slope'])
    print("ptemp.max,ftemp.max=",ptemp.max(),ftemp.max())
    
    #axscor.plot(ftemp,ptemp,'C0')
    
    
    #return

    ###########################################################################
    # LOOP OVER ROIS AND PLOT SCATTER IN APPROPRIATE COLOR
    ###########################################################################
    for R in ROI:
        RLatLims=-LatLims[0]+np.array([ROI[R][0],ROI[R][1]])
        RCM=ROI[R][2]
        RLonRng=ROI[R][3]
        RLonLims=[360-int(RCM+RLonRng),360-int(RCM-RLonRng)]

        RLonLims=np.array(np.array(RLonLims)-LonLims[0]).astype(int) #Just some extra protection

        print(RLatLims[0],RLatLims[1],
                         RLonLims[0],RLonLims[1])
        print(RLatLims[0]*scale,RLatLims[1]*scale,
                         RLonLims[0]*scale,RLonLims[1]*scale)
        
        subpatch1=patch1[RLatLims[0]*scale:RLatLims[1]*scale,
                         RLonLims[0]*scale:RLonLims[1]*scale]
        subpatch2=patch2[RLatLims[0]*scale:RLatLims[1]*scale,
                         RLonLims[0]*scale:RLonLims[1]*scale]*fNH3factor

        if dataversion=="H":
            axscor.scatter(subpatch2,subpatch1,marker=".",s=2,color=ROIcolors[R],linewidths=0,alpha=1.0,label=R)
        else:
            axscor.scatter(subpatch2,subpatch1,marker="o",s=3.0,color=ROIcolors[R],alpha=0.8,label=R)
            
        ROIout=statistics_helper(ROIout,obskey,R,subpatch1,subpatch2,
                                 ROIcolors[R],axscor,alpha=1.0,
                                 plot_ellipse=plot_ellipse)

    parent_results=mahalanobis_to_parent(ROIout, obskey)
    Mahal_pairwise=pairwise_mahalanobis(ROIout, obskey)
    Mahal_pairwise_roi=roi_pairwise_mahalanobis(ROIout, obskey)
    Mahalanobis_out={'parent_results':parent_results,
                     'pairwise':Mahal_pairwise,
                     'pairwise_roi':Mahal_pairwise_roi}   
    
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
    # 2. Grab what the automatic routine detected (scatter, plots, etc.)
    #handles, labels = axscor.get_legend_handles_labels()
    
    # 3. Check if the deep routine attached our custom ellipse proxy attribute
    #if hasattr(axscor, "ellipse_proxy"):
    #    handles.append(axscor.ellipse_proxy)
    #    labels.append(axscor.ellipse_proxy.get_label())
    
    # 4. Generate the complete, combined legend safely at the top level
    #axscor.legend(handles=handles, labels=labels,fontsize=8,ncols=3,labelcolor='mfc')                 
    axscor.legend(fontsize=7,ncols=4,labelcolor='mfc')
        
    return ROIout,Mahalanobis_out
  
