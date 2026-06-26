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
    widthi = 2 * np.sqrt(eigenvaluesi[1])
    heighti = 2 * np.sqrt(eigenvaluesi[0])
    # Angle of rotation in degrees (from the first eigenvector)
    anglei = np.degrees(np.arctan2(eigenvectorsi[1, 0], -eigenvectorsi[0, 0]))

    ellipse=Ellipse(
        xy=(mean_patch2, mean_patch1),
        width=widthi * 1.96,
        height=heighti * 1.96,
        angle=anglei,
        edgecolor=clr,
        alpha=alpha,
        facecolor="none",
        linestyle="--",
        linewidth=1.5)#,
        #label="95% Confidence")
    axscor.add_patch(ellipse)
    ## 2. Create a clean Line2D proxy matching your ellipse styling
    #axscor.ellipse_proxy = mlines.Line2D(
    #    [], [], color=clr, linestyle="--", linewidth=1.5, label="95% Confidence")   
    
def statistics_helper(ROIout,obskey,R,patch1,patch2,clr,axscor,alpha=1.0):
    mean_patch1,mean_patch2,stdv_patch1,stdv_patch2 = np.mean(patch1), np.mean(patch2),np.std(patch1),np.std(patch2)
    ROIout[obskey]['mean1'].append(mean_patch1)
    ROIout[obskey]['mean2'].append(mean_patch2)
    ROIout[obskey]['stdv1'].append(stdv_patch1)
    ROIout[obskey]['stdv2'].append(stdv_patch2)
    ROIout[obskey]['roilabel'].append(R)
    ROIout[obskey]['nsamples'].append(patch1.size)
    
    slopei, intercepti, r_valuei, p_valuei, std_erri = stats.linregress(patch1.ravel(), patch2.ravel())
    ROIout[obskey]['slope'].append(slopei)
    ROIout[obskey]['intercept'].append(intercepti)
    ROIout[obskey]['r_value'].append(r_valuei)
    ROIout[obskey]['p_value'].append(p_valuei)
    ROIout[obskey]['std_err'].append(std_erri)

    cov_matrixi = np.cov([patch1.ravel(), patch2.ravel()],rowvar=True)
    ROIout[obskey]['cov_matrix'].append(cov_matrixi)
    #print("cov_matrixi",cov_matrixi)
    plot_Mahal_ellipse(cov_matrixi,mean_patch1,mean_patch2,axscor,clr,alpha=0.8)
    
    return ROIout


    
def plot_roi_scatter(obskey,dateobs,ROI_ID,patch1,patch2,Real_CM2,LatLims,LonLims,axscor,PCldlow,PCldhigh,
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
    import numpy as np
    
    ###########################################################################
    # LOOP OVER ROIS AND PLOT SCATTER IN APPROPRIATE COLOR
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

    ROIout=statistics_helper(ROIout,obskey,'All',patch1,patch2,'grey',axscor,alpha=0.2)
    #print(np.array(ROIout[obskey]['cov_matrix'])[0,:,:])
    #plot_Mahal_ellipse(np.array(ROIout[obskey]['cov_matrix'])[0,:,:],
    #                   ROIout[obskey]['mean1'],
    #                   ROIout[obskey]['mean2'],axscor,'grey',alpha=0.8)

    #counter=1
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

        if dataversion=="H":
            axscor.scatter(subpatch2,subpatch1,marker=".",s=2,color=ROIcolors[R],linewidths=0,alpha=1.0,label=R)
        else:
            axscor.scatter(subpatch2,subpatch1,marker="o",s=3.0,color=ROIcolors[R],alpha=0.8,label=R)
            
        ROIout=statistics_helper(ROIout,obskey,R,subpatch1,subpatch2,ROIcolors[R],axscor,alpha=1.0)
        #print(counter, np.array(ROIout[obskey]['cov_matrix'])[counter,:,:])

        #plot_Mahal_ellipse(np.array(ROIout[obskey]['cov_matrix'])[counter,:,:],
        #                   ROIout[obskey]['mean1'][counter],
        #                   ROIout[obskey]['mean2'][counter],axscor,'grey',alpha=0.8)
        #counter+=1

    parent_results=mahalanobis_to_parent(ROIout, obskey)
    Mahal_pairwise=pairwise_mahalanobis(ROIout, obskey)
    Mahal_pairwise_roi=roi_pairwise_mahalanobis(ROIout, obskey)

    import json
    # Using a context manager ensures the file closes automatically
    GMM_file='C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/mahalanobis_clusters_with_covariances.json'
    with open(GMM_file, "r", encoding="utf-8") as f:
        GMM_Mahalanobis = json.load(f)
    
    total_clusters=6 #needs to be an input and so does the ROI
    if 'Ammonia' in xaxistitle:
        plottype='NH3_PCld'
        param1,param2='PCld','NH3'
    elif 'Color' in xaxistitle:
        plottype='AOI_CI'
        param1,param2='AOI','CI'
    else:
        plottype=False
    if plottype:
        subdict=GMM_Mahalanobis[obskey+ROI_ID][plottype][str(total_clusters)]
        for cluster_number in range(0,int(total_clusters)):
            mean1=float(subdict[param1][str(cluster_number+1)]['mean'])
            mean2=float(subdict[param2][str(cluster_number+1)]['mean'])
            print("GMM Data*************************************************")
            print(mean1,mean2)
            Test_Mahal_GMM_covariance=np.array(json.loads(subdict['covariances']))[cluster_number,:,:]
            print(Test_Mahal_GMM_covariance)
            print("*********************************************************")
            if plottype=='NH3_PCld':
                plot_Mahal_ellipse(np.flip(Test_Mahal_GMM_covariance),mean1,mean2,
                                   axscor,'C'+str(cluster_number),alpha=0.8)
            elif plottype=='AOI_CI':
                plot_Mahal_ellipse(Test_Mahal_GMM_covariance,mean1,mean2,
                                   axscor,'C'+str(cluster_number),alpha=0.8)
            axscor.scatter([],[],label='GMM '+str(cluster_number+1))
        """
    print()
    print("############## Mahalanobis to parent")
    for roi, vals in parent_results.items():
        print(
            roi,
            vals['D2'],
            vals['D']
        )
    print("############## Mahalanobis pairwise, pooled")
    print()
    
    print(Mahal_pairwise['labels'],Mahal_pairwise['D'])
    print()
    print("############## Mahalanobis pairwise, pooled 4")
    print()
    print(Mahal_pairwise_roi['labels'],Mahal_pairwise_roi['D'])
    print("##############")
    print()
    """
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
    Mahalanobis_out={'parent_results':parent_results,
                     'pairwise':Mahal_pairwise,
                     'pairwise_roi':Mahal_pairwise_roi}   
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
  
