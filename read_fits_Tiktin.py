from config_VA import config_VA

def read_fits_L3_V2_helper(pathandfile,target="Jupiter",LonSys='3',
                     LimbCorrection=False):
    
    import get_spice_ephem as sp_ephem
    from astropy.io import fits
    import numpy as np
    from copy import deepcopy

    dataobj={}
    
    hdulist=fits.open(pathandfile)
    hdulist.info()
    dataobj['hdr']=hdulist[0].header
    dataobj['data']=np.flipud(hdulist[0].data[0])
    dataobj['lon']=np.flipud(hdulist[1].data)
    dataobj['lat']=np.flipud(hdulist[2].data)
    dataobj['sza']=np.flipud(hdulist[3].data) #!!!! Might have eza and sza reversed
    dataobj['eza']=np.flipud(hdulist[4].data) #!!!! Might have eza and sza reversed
    hdulist.close()

    CM3ck=float(dataobj['hdr']["HIERARCH PLANMAP SUBPOINT LON"])
    eph=sp_ephem.get_spice_ephem(dataobj['hdr']["DATE-OBS"],planet='target')
    dataobj['CM1']=float(eph[0].strip())
    dataobj['CM2']=float(eph[1].strip())
    dataobj['CM3']=float(eph[2].strip())
    if np.abs(CM3ck-dataobj['CM3'])>0.01:
        print("CM ERROR")
        
    if LimbCorrection and 'BUNIT' in dataobj['hdr']:
        exponent={'Cloud-top Press':0.25,'Mole Fraction':0.55}
        amfdata=(1.0/(np.cos(dataobj['sza']*np.pi/180))+1.0/np.cos(dataobj['eza']*np.pi/180.))/2.0
        dataobj['data']=dataobj['data']*(amfdata**exponent[dataobj['hdr']['BUNIT']])
        
    dataobj['datar']=deepcopy(np.roll(dataobj['data'],int(dataobj['CM3']-dataobj['CM'+LonSys]),axis=1))
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
                        LimbCorrection=False,pathin=config_VA[2]):

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

    #pathSci=pathin+"New_Results/"+obskey[:-1]+"/"+obskey+"/"+Level+'/'
    pathSci=pathin+"/"+obskey[:-1]+"/"+obskey+"/"+Level+'/'
    filesSci=os.listdir(pathSci)
    #pathIGB=pathin+'New_Results/'+obskey[:-1]+'/'+obskey+'/L1/'
    pathIGB=pathin+'/'+obskey[:-1]+'/'+obskey+'/L1/'
    filesIGB=os.listdir(pathIGB)
    
    sciobjects={'PCld':[],'fNH3':[]}
    RGBobjects={'NIR':[],'GRN':[],'BLU':[]}
    for key in sciobjects:
        filename=[item for item in filesSci if key in item][0]
        sciobjects[key]=read_fits_L3_V2_helper(pathSci+filename,target="Jupiter",
                             LonSys=LonSys,LimbCorrection=LimbCorrection)
    
    for key in RGBobjects:
        print(key)
        filename=[item for item in filesIGB if key in item][0]
        RGBobjects[key]=read_fits_L3_V2_helper(pathIGB+filename,target="Jupiter",
                             LonSys=LonSys,LimbCorrection=LimbCorrection)

    IGBdatar=np.dstack((RGBobjects['NIR']['datar'],RGBobjects['GRN']['datar'],RGBobjects['BLU']['datar']))
    IGBdatarx=np.nan_to_num(IGBdatar, nan=0.0, posinf=1.0, neginf=0.0)

    return(sciobjects['PCld']['hdr'],sciobjects['PCld']['datar'],sciobjects['PCld']['CM'+LonSys],
           sciobjects['fNH3']['hdr'],sciobjects['fNH3']['datar'],sciobjects['fNH3']['CM'+LonSys],
           IGBdatarx/np.max(IGBdatarx),RGBobjects['GRN']['CM'+LonSys],RGBobjects['GRN']['hdr']['DATE-OBS'])

def read_fits_map_L3_V1(obskey="20231026UTa",imagetype="Map",Level="L3",
                        target="Jupiter",LonSys='3',
                        LimbCorrection=False,pathin=config_VA[1]):
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
