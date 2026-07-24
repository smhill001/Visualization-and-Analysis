def flatten_patch(patch,limitx=True,limity=True):
    """
    Empirical longitudinal flattening of a map patch using a 2nd order
    polynomial fit

    Parameters
    ----------
    patch : TYPE
        DESCRIPTION.

    Returns
    -------
    patchflat : TYPE
        DESCRIPTION.

    """
    import numpy as np
    x=np.arange(0,patch.shape[1],1)
    if limitx:
        x0=int(patch.shape[1]/6.)
        x1=int(patch.shape[1]*5./6.)
        xlim=x[x0:x1]
    else:
        x0=0
        x1=patch.shape[1]
        xlim=x
    if limity:
        y0=int(patch.shape[0]/6.)
        y1=int(patch.shape[0]*5./6.)
        #ylim=x[y0:y1]
    else:
        y0=0
        y1=patch.shape[0]
        #ylim=y
    
    
    
    #y=np.mean(patch,axis=0) #Meridional average
    ylim=np.mean(patch[y0:y1,x0:x1],axis=0) #Meridional average
    
    coefs=np.polyfit(xlim,ylim,2)
    linfit=coefs[2]+coefs[1]*x+coefs[0]*x**2
    
    linfitnorm=linfit/np.mean(linfit)
    
    patchflat=patch/linfitnorm
    
    return patchflat