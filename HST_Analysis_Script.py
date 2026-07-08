# -*- coding: utf-8 -*-

def HST_Analysis_Script(obskeyHST,LonSys,HST=True,SCT=False,ROI=False,ROI_ID='',
                        compare=False,segment=False,fNH3factor=1.0):
    
    import sys
    import socket
    hostname = socket.gethostname()
    from config_VA import Host_path
    from HST_Study_Maps import Study_Maps
    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/HST/')
    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')
    #sys.path.append()
    import L3_Jup_Map_Plot_V2 as L3MP

    collection=Study_Maps[obskeyHST]['collection']
    plotoptions=Study_Maps[obskeyHST][LonSys]['plotoptions']
    if ROI_ID in Study_Maps[obskeyHST][LonSys]:
        ROI=Study_Maps[obskeyHST][LonSys][ROI_ID]['ROI']
        CoLatLims=Study_Maps[obskeyHST][LonSys][ROI_ID]['CoLatLims']
        LonRng=Study_Maps[obskeyHST][LonSys][ROI_ID]['LonRng']
        CMpref=Study_Maps[obskeyHST][LonSys][ROI_ID]['CMpref']
    else:
        CoLatLims=Study_Maps[obskeyHST][LonSys]['CoLatLims']
        LonRng=Study_Maps[obskeyHST][LonSys]['LonRng']
        CMpref=Study_Maps[obskeyHST][LonSys]['CMpref']
            
    if HST:
        print("@@@@@@@@@@@@@@@@@@@@@ ROI_ID=",ROI_ID)
        L3MP.L3_Jup_Map_Plot_V2(obskey=obskeyHST, 
                                CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                                subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                                dataversion='H',smoothcont=5,segment=segment,ROI_ID=ROI_ID,
                                ROI=ROI,compare=compare,fNH3factor=fNH3factor)

    if SCT:
        L3MP.L3_Jup_Map_Plot_V2(obskey=collection, 
                                CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                                subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                                dataversion=2,smoothcont=0,segment=segment,ROI_ID=ROI_ID,
                                ROI=ROI,compare=compare,fNH3factor=fNH3factor)
