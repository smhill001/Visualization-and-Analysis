import numpy as np
def make_lat_lon_str(lats,LonLimsEast):
    """
    Creates strings for insertion into file names for the latitude and
    longitude range of a give map.

    Parameters
    ----------
    lats : numpy array
        DESCRIPTION.
    LonLims : numpy array
        DESCRIPTION.

    Returns
    -------
    latstr : string
        Latitude range description string for use in filenames.
    lonstr : string
        Longitude range description string for use in filenames.

    """
    if int(lats[0])<90:
        latstr="N"+str(90-lats[0])
    if int(lats[0])==90:
        latstr="S0"+str(90-lats[0])
    if int(lats[0])>90:
        latstr="S"+str(lats[0]-90)
        
    if int(lats[1])<90:
        latstr=latstr+"-N"+str(90-lats[1])
    if int(lats[1])==90:
        latstr=latstr+"-S0"+str(90-lats[1])
    if int(lats[1])>90:
        latstr=latstr+"-S"+str(lats[1]-90)
        
    lonstr=str(np.mod(360-LonLimsEast[1],360)).zfill(3)+"-"+str(np.mod(360-LonLimsEast[0],360)).zfill(3)
    
    return latstr,lonstr