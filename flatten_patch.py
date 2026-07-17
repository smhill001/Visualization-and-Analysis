def flatten_patch(patch,limit=False):
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
    if limit:
        y0=int(patch.shape[0]/6.)
        y1=int(patch.shape[0]*5./6.)
    else:
        y0=0
        y1=patch.shape[0]
    
    
    
    #y=np.mean(patch,axis=0) #Meridional average
    y=np.mean(patch[y0:y1,:],axis=0) #Meridional average
    
    coefs=np.polyfit(x,y,2)
    linfit=coefs[2]+coefs[1]*x+coefs[0]*x**2
    
    linfitnorm=linfit/np.mean(linfit)
    
    patchflat=patch/linfitnorm
    
    return patchflat