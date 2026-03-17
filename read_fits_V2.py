from config_VA import config_VA

def read_fits_L3_V2_helper(pathandfile,target="Jupiter",LonSys='3',
                     LimbCorrection=False,dataversion=2):
    
    import get_spice_ephem as sp_ephem
    from astropy.io import fits
    import numpy as np
    from copy import deepcopy

    dataobj={}
    
    hdulist=fits.open(pathandfile)
    hdulist.info()
    dataobj['hdr']=hdulist[0].header
    if dataversion=='H':
        #dataobj['data']=np.flipud(hdulist[0].data)
        dataobj['data']=hdulist[0].data
    else:
        dataobj['data']=np.flipud(hdulist[0].data[0])
        
    if dataversion!='H':
        dataobj['lon']=np.flipud(hdulist[1].data)
        dataobj['lat']=np.flipud(hdulist[2].data)
        dataobj['sza']=np.flipud(hdulist[3].data) #!!!! Might have eza and sza reversed
        dataobj['eza']=np.flipud(hdulist[4].data) #!!!! Might have eza and sza reversed
    hdulist.close()

    eph=sp_ephem.get_spice_ephem(dataobj['hdr']["DATE-OBS"],planet='target')
    dataobj['hdr']['CM1']=float(eph[0].strip())
    dataobj['hdr']['CM2']=float(eph[1].strip())
    dataobj['hdr']['CM3']=float(eph[2].strip())
    print(dataobj['hdr']['CM1'])
    print(dataobj['hdr']['CM2'])
    print(dataobj['hdr']['CM3'])
    #CM3ck=float(dataobj['hdr']["HIERARCH PLANMAP SUBPOINT LON"])
    #if np.abs(CM3ck-dataobj['hdr']['CM3'])>0.01:
    #    print("CM ERROR")
        
    if LimbCorrection and 'BUNIT' in dataobj['hdr']:
        exponent={'Cloud-top Press':0.25,'Mole Fraction':0.55}
        amfdata=(1.0/(np.cos(dataobj['sza']*np.pi/180))+1.0/np.cos(dataobj['eza']*np.pi/180.))/2.0
        dataobj['data']=dataobj['data']*(amfdata**exponent[dataobj['hdr']['BUNIT']])
        
    print(dataobj['data'].shape)
    if dataversion==1 or dataversion==2:
        dataobj['datar']=deepcopy(np.roll(dataobj['data'],int(dataobj['hdr']['CM3']-dataobj['hdr']['CM'+LonSys]),axis=1))
    elif dataversion=='H':
        #dataobj['datar']=deepcopy(np.roll(dataobj['data'],20*int(dataobj['hdr']['CM3']-dataobj['hdr']['CM'+LonSys]),axis=1))
        #!!!! Below is for customization requiring that HST data be read from a FITS file of the same LonSys as will be plotted
        dataobj['datar']=deepcopy(np.roll(dataobj['data'],20*int(dataobj['hdr']['CM'+LonSys]-dataobj['hdr']['CM'+LonSys]),axis=1))
    """
    fig1,axs1=pl.subplots(1,figsize=(3,6), dpi=150, facecolor="white")
    fig2,axs2=pl.subplots(1,figsize=(3,6), dpi=150, facecolor="white")
    print("############",dataobj['data'].shape,dataobj['datar'].shape,int(dataobj['CM3']-dataobj['CM'+LonSys]))
    axs1.imshow(dataobj['data'])
    axs2.imshow(dataobj['datar'])
    """
    return dataobj

def read_fits_map_L3_V2(obskey="20251016UTa",imagetype="Map",Level="L3",
                        target="Jupiter",LonSys='3',
                        LimbCorrection=False,dataversion=2):

    import sys
    import os
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Maps')
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Services')

    sys.path.append('./Services')
    import numpy as np
    pathin=config_VA[dataversion]
    #pathSci=pathin+"New_Results/"+obskey[:-1]+"/"+obskey+"/"+Level+'/'
    pathSci=pathin+"/"+obskey[:-1]+"/"+obskey+"/"+Level+'/'
    filesSciTemp=os.listdir(pathSci)
    print("############## filesSciTemp= ",filesSciTemp)

    filesSci=[]
    print("######## dataversion=",dataversion)
    if dataversion=='H':
        for fn in filesSciTemp:
            if 'Sys'+LonSys in fn:
                filesSci.append(fn)
    else:
        filesSci=filesSciTemp
    #pathIGB=pathin+'New_Results/'+obskey[:-1]+'/'+obskey+'/L1/'
    pathIGB=pathin+'/'+obskey[:-1]+'/'+obskey+'/L1/'
    filesIGBTemp=os.listdir(pathIGB)
    filesIGB=[]
    if dataversion=='H':
        for fn in filesIGBTemp:
            if 'Sys'+LonSys in fn:
                filesIGB.append(fn)
    else:
        filesIGB=filesIGBTemp

    sciobjects={'PCld':[],'fNH3':[]}
    if dataversion==2:
        RGBobjects={'NIR':[],'GRN':[],'BLU':[]}
    elif dataversion=='H':
        RGBobjects={'673':[],'502':[],'395':[]}

    print("############## filesSci= ",filesSci)

    for key in sciobjects: 
        filename=[item for item in filesSci if key in item][0]
        print(filename)
        sciobjects[key]=read_fits_L3_V2_helper(pathSci+filename,target="Jupiter",
                             LonSys=LonSys,LimbCorrection=LimbCorrection,dataversion=dataversion)
    
    for key in RGBobjects:
        print(key)    
        templist=[item for item in filesIGB if key in item]
        if not templist:
            filename=[item for item in filesIGB if RGBobjects[0] in item][0] #Kludge if missing GRN or BLU
            ##!!! Should use 'try' logic to trap error and find any working RGB channel
            RGBobjects[key]=read_fits_L3_V2_helper(pathIGB+filename,target="Jupiter",
                                 LonSys=LonSys,LimbCorrection=LimbCorrection,dataversion=dataversion)
        elif templist:
            filename=templist[0]
            RGBobjects[key]=read_fits_L3_V2_helper(pathIGB+filename,target="Jupiter",
                                 LonSys=LonSys,LimbCorrection=LimbCorrection,dataversion=dataversion)

    if dataversion==2:
        IGBdatar=np.dstack((RGBobjects['NIR']['datar'],RGBobjects['GRN']['datar'],RGBobjects['BLU']['datar']))
    elif dataversion=='H':
        IGBdatar=np.dstack((RGBobjects['673']['datar'],RGBobjects['502']['datar'],RGBobjects['395']['datar']))
    IGBdatarx=np.nan_to_num(IGBdatar, nan=0.0, posinf=1.0, neginf=0.0)

    
    if dataversion==2:
        return(sciobjects['PCld']['hdr'],sciobjects['PCld']['datar'],
               sciobjects['fNH3']['hdr'],sciobjects['fNH3']['datar'],
               IGBdatarx/np.max(IGBdatarx),RGBobjects['GRN']['hdr']['CM'+LonSys],RGBobjects['GRN']['hdr']['DATE-OBS'])
    elif dataversion=='H':
        return(sciobjects['PCld']['hdr'],sciobjects['PCld']['datar'],
               sciobjects['fNH3']['hdr'],sciobjects['fNH3']['datar'],
               IGBdatarx/np.max(IGBdatarx),RGBobjects['502']['hdr']['CM'+LonSys],RGBobjects['502']['hdr']['DATE-OBS'])

def read_fits_map_L3_V1(obskey="20231026UTa",imagetype="Map",Level="L3",
                        target="Jupiter",LonSys='3',
                        LimbCorrection=False,dataversion=1):
    """
    Created on Mon Nov 20 08:42:28 2023
    Called by: Map_Jup_Atm_2022_P3, currently only for L3 data to plot
               maps of fNH3 and PCld
               
    Updated 7/6/2025 by SMH to incorporate empirical limb correction (just for L3 for now)
    
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    from matplotlib.pyplot import imread
    from astropy.io import fits
    sys.path.append('./Services')
    import get_obs_list as getlist

    import get_spice_ephem as sp_ephem
    import numpy as np

    pathin=config_VA[dataversion]
    sourcedata=obskey#+"_"+imagetype
    sourcefiles=getlist.get_obs_list(planet=target)
    pathRGB='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obskey[0:10]+'/'
    pathmicron="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/FITS/output/"
          
    if imagetype=="Map":
        RGBfile=sourcefiles[sourcedata]["RGBfile"]+"_CM2_L360_MAP-BARE.png"
    elif imagetype=="Img":
        RGBfile=sourcefiles[sourcedata]["RGBfile"]+".png"

    #Only read an RGB file if there's one in the sourcefiles and obskey    
    if RGBfile != 'NA': 
        RGB=imread(pathRGB+RGBfile)
        RGBsec=str(int(str(RGBfile[16:17]))*6)
        RGBtime=(RGBfile[0:10]+"_"+RGBfile[11:13]+":"+RGBfile[13:15]+":"+RGBsec.zfill(2))
        eph=sp_ephem.get_spice_ephem(RGBtime,planet=target)
        RGB_CM1=float(eph[0].strip())
        RGB_CM2=float(eph[1].strip())
        RGB_CM3=float(eph[2].strip())

    #Create suffixes and paths for L3 FITS files of Ammonia and Clouds
    CH4suffix="-"+target+"_"+imagetype+"_L3PCld_S0"
    NH3suffix="-"+target+"_"+imagetype+"_L3fNH3_S0"
    pathFITS=pathin
    print(pathFITS)

    # try-except logic for when variation files existed. Hasn't been run 
    #   probably since early 2023 at the latest. Doubtful it would work now 
    #   (2024-11-18)
    try:
        PCloudfile=sourcefiles[sourcedata]['CH4file'][0:17]+CH4suffix+\
                sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Metadata']['Variation']
    except:
        PCloudfile=sourcefiles[sourcedata]['CH4file'][0:17]+CH4suffix+".fits"
        variation=""
    try:
        fNH3file=sourcefiles[sourcedata]['NH3file'][0:17]+NH3suffix+\
                sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Metadata']['Variation']
    except:
        fNH3file=sourcefiles[sourcedata]['NH3file'][0:17]+NH3suffix+".fits"
        variation=""
    
    #Read the files!!!
    PCloudhdulist=fits.open(pathFITS+PCloudfile)
    PCloudhdulist.info()
    PCloudhdr=PCloudhdulist[0].header
    PClouddata=PCloudhdulist[0].data
    PCloudhdulist.close()
    
    fNH3hdulist=fits.open(pathFITS+fNH3file)
    fNH3hdulist.info()
    fNH3hdr=fNH3hdulist[0].header
    fNH3data=fNH3hdulist[0].data
    fNH3sza=fNH3hdulist[1].data
    fNH3eza=fNH3hdulist[2].data
    fNH3hdulist.close()
    
    NH3time=fNH3hdr["DATE-OBS"]
    NH3_CM1=float(fNH3hdr["CM1"])
    NH3_CM2=float(fNH3hdr["CM2"])
    NH3_CM3=float(fNH3hdr["CM3"])
        
    CH4time=PCloudhdr["DATE-OBS"]
    CH4_CM1=float(PCloudhdr["CM1"])
    CH4_CM2=float(PCloudhdr["CM2"])
    CH4_CM3=float(PCloudhdr["CM3"])
    print("***********LonSys=",LonSys)

    if LonSys=='1':
        #print("LonSys1=",LonSys)
        if imagetype=="Map":
            RGBroll=RGB_CM2-RGB_CM1
            RGB=np.roll(RGB,int(RGBroll),axis=1)

        NH3roll=NH3_CM3-NH3_CM1
        fNH3datar=np.roll(fNH3data,int(NH3roll),axis=1)
        
        CH4roll=CH4_CM3-CH4_CM1
        #print("CH4roll=",CH4roll)
        PClouddatar=np.roll(PClouddata,int(CH4roll),axis=1)
        
        fNH3szar=np.roll(fNH3sza,int(NH3roll),axis=1)
        fNH3ezar=np.roll(fNH3eza,int(NH3roll),axis=1)
        
        NH3CM=NH3_CM1
        CH4CM=CH4_CM1
        RGBCM=RGB_CM1       
        #CM=NH3_CM1
        #Real_CM=NH3_CM1

    if LonSys=='2':
        #print("In LonSys==2")
        NH3roll=NH3_CM3-NH3_CM2
        fNH3datar=np.roll(fNH3data,int(NH3roll),axis=1)
        
        CH4roll=CH4_CM3-CH4_CM2
        #print("&&&&&&&&&& CH4roll=",CH4roll)
        PClouddatar=np.roll(PClouddata,int(CH4roll),axis=1)

        fNH3szar=np.roll(fNH3sza,int(NH3roll),axis=1)
        fNH3ezar=np.roll(fNH3eza,int(NH3roll),axis=1)
        
        NH3CM=NH3_CM2
        CH4CM=CH4_CM2
        RGBCM=RGB_CM2
        
    if LonSys=='3':
        if imagetype=="Map":
            RGBroll=RGB_CM3-RGB_CM2
            RGB=np.roll(RGB,int(-RGBroll),axis=1)

        fNH3datar=fNH3data
        PClouddatar=PClouddata
        fNH3szar=fNH3sza
        fNH3ezar=fNH3eza

        NH3CM=NH3_CM3
        CH4CM=CH4_CM3
        RGBCM=RGB_CM3
            
    if LimbCorrection:
        amfdata=(1.0/fNH3szar+1.0/fNH3ezar)/2.0
        fNH3datar=fNH3datar*(amfdata**0.55)
        PClouddatar=PClouddatar*amfdata**0.25
    
    return(PCloudhdr,PClouddatar,
           fNH3hdr,fNH3datar,
           fNH3szar,fNH3ezar,
           RGB,RGBCM,RGBtime)

def read_fits_map_L4(LonSys,collection="20220904-20220905",
                      collectionIRTF="20240205-20240205",
                      collection889CH4="20240131-20240201"):
    """
    Retrieve continguous maps from FITS files for plotting and analysis.
    
    !!!Development notes: 
    1) Currently assumes that fNH3, PCld, and RGB 'optical' data is always
    provided. It could be useful to make that optional like the IRTF and 889-nm
    data in order to generalize this function.
    2) A helper function to read the m x n FITS files and perhaps the RGB file
    could dramatically shorten this function and make it more readable.
    

    Parameters
    ----------
    LonSys : String
        Indicates desired longitude system of map plot
    collection : String, optional
        The map collection label for contiguous map FITS files of fNH3, PCld,
        and visual context
    collectionIRTF : String, optional
        The map collection label for contiguous map FITS file of 5-micron data
    collection889CH4 : String, optional
        The map collection label for contiguous map FITS file of 889-nm data

    Returns
    -------
    fNH3data : float array
        Cylindrically mapped mean ammonia mole fraction at each lat-lon in ppm.
    fNH3stdv : float array
        Cylindrically mapped fNH3 standard deviation at each lat-lon in ppm.
    fNH3frac : float array
        Cylindrically mapped fractional uncertainty in fNH3, e.g., 
        fNH3stdv/fNH3data
    fNH3time : float array
        Cylindrically mapped mean fNH3 observation time at each lat-lon in JD(?).
    PClddata : float array
        Cylindrically mapped mean effective cloud pressure at each lat-lon in 
        mbar.
    PCldstdv : float array
        Cylindrically mapped PCld standard deviation at each lat-lon in mbar.
    PCldfrac : float array
        Cylindrically mapped fractional uncertainty in PCld, e.g., 
        PCldstdv/PClddata
    PCldtime : float array
        Cylindrically mapped mean PCld observation time at each lat-lon in JD(?).
    RGBdata : TYPE
        DESCRIPTION.
    RGBstdv : TYPE
        DESCRIPTION.
    RGBfrac : TYPE
        DESCRIPTION.
    RGBtime : TYPE
        DESCRIPTION.
    IRTFdata : float array
        Cylindrically mapped mean effective cloud pressure at each lat-lon in 
        uncalibrated units.
    IRTFstdv : float array
        Cylindrically mapped IRTF standard deviation at each lat-lon in mbar.
    IRTFfrac : float array
        Cylindrically mapped fractional uncertainty in IRTF, e.g., 
        IRTFstdv/IRTFdata
    IRTFtime : float array
        Cylindrically mapped mean IRTF observation time at each lat-lon in JD(?).
    JALPOdata : TYPE
        DESCRIPTION.

    """
    from astropy.io import fits
    import os
    import numpy as np
    print()
    print("############ IN READ Contiguous MAP")
    print("####### LonSys,collection,IRTFcollection,CH4889collection",
          LonSys,collection,collectionIRTF,collection889CH4)
    print()

    
    cm_key = {
        '1': 'CM1',
        '2': 'CM2',
        '3': 'CM3'#,
        #'subobs':180
    }.get(LonSys)
    
    ctype1_key = {
        'Sys. 1 Longitude':'CM1',
        'Sys. 2 Longitude':'CM2',
        'Sys. 3 Longitude':'CM3',
        }
    
    ###########################################################################
    # Identify files to read
    ###########################################################################
    pathinp="C:/Astronomy/Projects/SAS 2021 Ammonia/Data/L4 FITS (cont maps)/"
    contents = os.listdir(pathinp)
    files_in_directory = [item for item in contents if os.path.isfile(os.path.join(pathinp, item))]
    
    fits_in_directory = [item for item in files_in_directory \
                         if '.fits' in item ]  
    collection_in_directory = [item for item in fits_in_directory \
                         if collection in item]  
    LonSys_in_directory = [item for item in collection_in_directory \
                         if 'Sys'+LonSys in item]
    print("############### LonSys_in_directory=",LonSys_in_directory)

    ###########################################################################
    # Read fNH3 data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    fNH3file = [item for item in LonSys_in_directory \
                         if 'L4fNH3' in item]  
    print("############### fNH3file=",fNH3file)
    if len(fNH3file)>0:
        print("############### fNH3file=",fNH3file[0])
        fNH3hdulist=fits.open(pathinp+fNH3file[0])
        fNH3hdulist.info()
        fNH3hdr=fNH3hdulist[0].header
        roll=int(fNH3hdr[ctype1_key[fNH3hdr['CTYPE1']]]-fNH3hdr[cm_key])
        fNH3data=np.roll(fNH3hdulist[0].data,roll,axis=1)
        fNH3stdv=np.roll(fNH3hdulist[1].data,roll,axis=1)
        fNH3frac=np.roll(fNH3hdulist[2].data,roll,axis=1)
        fNH3time=np.roll(fNH3hdulist[3].data,roll,axis=1)
        fNH3hdulist.close()       
        
    ###########################################################################
    # Read PCld data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    PCldfile = [item for item in LonSys_in_directory \
                         if 'L4PCld' in item]      
    if len(PCldfile)>0:
        PCldhdulist=fits.open(pathinp+PCldfile[0])
        PCldhdulist.info()
        PCldhdr=PCldhdulist[0].header
        roll=int(PCldhdr[ctype1_key[PCldhdr['CTYPE1']]]-PCldhdr[cm_key])
        PClddata=np.roll(PCldhdulist[0].data,roll,axis=1)
        PCldstdv=np.roll(PCldhdulist[1].data,roll,axis=1)
        PCldfrac=np.roll(PCldhdulist[2].data,roll,axis=1)
        PCldtime=np.roll(PCldhdulist[3].data,roll,axis=1)
        PCldhdulist.close()
        
    ###########################################################################
    # Read RGB data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    RGBfile = [item for item in LonSys_in_directory \
                         if 'L4RGB' in item]      
    if len(RGBfile)>0:
        RGBhdulist=fits.open(pathinp+RGBfile[0])
        RGBhdulist.info()
        RGBhdr=RGBhdulist[0].header
        roll=int(RGBhdr[ctype1_key[RGBhdr['CTYPE1']]]-RGBhdr[cm_key])
        RGBdata=np.roll(RGBhdulist[0].data,roll,axis=1)
        RGBstdv=np.roll(RGBhdulist[1].data,roll,axis=1)
        RGBfrac=np.roll(RGBhdulist[2].data,roll,axis=1)
        RGBtime=np.roll(RGBhdulist[3].data,roll,axis=1)
        RGBhdulist.close()

    ###########################################################################
    # Read IRTF data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    if collectionIRTF:        
        collectionIRTF_in_directory = [item for item in fits_in_directory \
                             if collectionIRTF in item]  
        IRTFfile = [item for item in collectionIRTF_in_directory \
                             if 'L4IRTF' in item]      
        if len(IRTFfile)>0:
            print("############### IRTF Map=",IRTFfile)
            IRTFhdulist=fits.open(pathinp+IRTFfile[0])
            IRTFhdulist.info()
            IRTFhdr=IRTFhdulist[0].header
            roll=int(IRTFhdr['CM3']-IRTFhdr[cm_key])
            IRTFdata=np.roll(IRTFhdulist[0].data,roll,axis=1)
            IRTFstdv=np.roll(IRTFhdulist[1].data,roll,axis=1)
            IRTFfrac=np.roll(IRTFhdulist[2].data,roll,axis=1)
            IRTFtime=np.roll(IRTFhdulist[3].data,roll,axis=1)
            IRTFhdulist.close()
        else:
            IRTFdata=np.zeros((180,360))
            IRTFstdv=np.zeros((180,360))
            IRTFfrac=np.zeros((180,360))
            IRTFtime=np.zeros((180,360))       
    else:
        IRTFdata=np.zeros((180,360))
        IRTFstdv=np.zeros((180,360))
        IRTFfrac=np.zeros((180,360))
        IRTFtime=np.zeros((180,360))

    ###########################################################################
    # Read 889-nm JALPO data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    if collection889CH4:
        
        collection889CH4_in_directory = [item for item in fits_in_directory \
                             if collection889CH4 in item]  
        JALPOfile = [item for item in collection889CH4_in_directory \
                             if 'L4889CH4' in item]      
        if len(JALPOfile)>0:
            JALPOhdulist=fits.open(pathinp+JALPOfile[0])
            JALPOhdulist.info()
            JALPOhdr=JALPOhdulist[0].header
            print("########### ",JALPOhdr['DATE-OBS'])
            print("########### ",JALPOhdr['CM3'])
            print("########### ",cm_key)
    
            roll=int(JALPOhdr['CM3'])-int(JALPOhdr[cm_key])
            JALPOdata=np.roll(JALPOhdulist[0].data,roll,axis=1)
            JALPOhdulist.close()
        else:
            JALPOdata=np.zeros((180,360))
    else:
        JALPOdata=np.zeros((180,360))
    
    return fNH3data,fNH3stdv,fNH3frac,fNH3time,fNH3hdr, \
        PClddata,PCldstdv,PCldfrac,PCldtime,PCldhdr, \
        RGBdata,RGBstdv,RGBfrac,RGBtime,RGBhdr, \
        IRTFdata,IRTFstdv,IRTFfrac,IRTFtime, \
        JALPOdata