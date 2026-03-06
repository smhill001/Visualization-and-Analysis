def flatten_patch(patch):
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
    y=np.mean(patch,axis=0)
    
    coefs=np.polyfit(x,y,2)
    linfit=coefs[2]+coefs[1]*x+coefs[0]*x**2
    
    linfitnorm=linfit/np.mean(linfit)
    
    patchflat=patch/linfitnorm
    
    return patchflat