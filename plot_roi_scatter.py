import numpy as np
from scipy import stats
from matplotlib.patches import Ellipse

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

    return labels, D

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

    for i in range(n):
        for j in range(i + 1, n):

            d = means[i] - means[j]

            D2 = d.T @ Sinv @ d

            D[i, j] = np.sqrt(D2)
            D[j, i] = D[i, j]

    return labels, D

def statistics_helper(ROIout,obskey,R,patch1,patch2,clr,axscor,alpha=1.0):
    ROIout[obskey]['mean1'].append(np.mean(patch1))
    ROIout[obskey]['mean2'].append(np.mean(patch2))
    ROIout[obskey]['stdv1'].append(np.std(patch1))
    ROIout[obskey]['stdv2'].append(np.std(patch2))
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
        edgecolor=clr,
        alpha=alpha,
        facecolor="none",
        linestyle="--",
        linewidth=1.5,
        label="95% confidence")
    axscor.add_patch(ellipse)

    return ROIout      
    
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
        
    parent_results=mahalanobis_to_parent(ROIout, obskey)
    
    print()
    print("##############")
    for roi, vals in parent_results.items():
        print(
            roi,
            vals['D2'],
            vals['D']
        )
    print("##############")
    print()
   
    labels,D=pairwise_mahalanobis(ROIout, obskey)
    print(labels,D)
    print("##############")
    print()
    
    labels4,D4=roi_pairwise_mahalanobis(ROIout, obskey)
    print(labels4,D4)
    print("##############")
    print()

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
        
    return ROIout
    #return(roilabel,mean1,stdv1,mean2,stdv2,slope,intercept,r_value,p_value,std_err,cov_matrix)#,meanamf)
  
