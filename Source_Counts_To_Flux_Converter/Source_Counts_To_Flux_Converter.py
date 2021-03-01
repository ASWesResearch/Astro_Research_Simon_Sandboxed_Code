import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
from astropy.io.fits import Header
from astropy.io import fits
import os
from astroquery.ned import Ned
import astropy.units as u
import sys
from astropy.io import ascii
def Evt2_File_Query(ObsID):
    query_path='/Volumes/xray/simon/all_chandra_observations/'+str(ObsID)+'/primary/*evt2*'
    evtfpath_L=glob.glob(query_path)
    if(len(evtfpath_L)!=1):
        #print str(ObsID)+" has "+str(len(evtfpath_L)-1)+"Additional Evt2 Files ! ! !"
        print str(ObsID)+" has "+str(len(evtfpath_L))+" Evt2 Files ! ! !"
        return "Error"
    evtfpath=evtfpath_L[0]
    return evtfpath
def Gname_Query(ObsID,NED_Resolved_Bool=False):
    #query_path='/Volumes/xray/simon/all_chandra_observations/'+str(ObsID)+'/primary/*evt2*'
    if(NED_Resolved_Bool==False):
        Flux_90_Fpath_L=glob.glob("/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/*/Flux_90_Files/"+str(ObsID)+"/*_Flux_90.txt")
        if(len(Flux_90_Fpath_L)!=1):
            #print str(ObsID)+" has "+str(len(Flux_90_Fpath_L)-1)+"Additional Flux 90 Files ! ! !"
            print str(ObsID)+" has "+str(len(Flux_90_Fpath_L))+" Flux 90 Files ! ! !"
            return "Error"
        Flux_90_Fpath=Flux_90_Fpath_L[0]
        Gname=Flux_90_Fpath.split("/")[7]
        #print "Gname: ", Gname
        return Gname
    if(NED_Resolved_Bool):
        Obs_ID=int(ObsID)
        path_Source_Flux_Table=os.path.realpath('/Volumes/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_Table.csv') #Absolute Path
        #print "path_Source_Flux_Table : ", path_Source_Flux_Table
        Source_Flux_Table=ascii.read(path_Source_Flux_Table)
        #print "Source_Flux_Table : \n", Source_Flux_Table
        Obs_ID_A=Source_Flux_Table["OBSID"]
        Obs_ID_L=list(Obs_ID_A)
        #print "Obs_ID_L : ", Obs_ID_L
        Obs_ID_Inx=Obs_ID_L.index(Obs_ID)
        #print "Obs_ID_Inx : ", Obs_ID_Inx
        Resolved_Name_A=Source_Flux_Table["resolvedObject"]
        #print "Resolved_Name_A : ", Resolved_Name_A
        Resolved_Name_L=list(Resolved_Name_A)
        #print "Resolved_Name_L : ", Resolved_Name_L
        Resolved_Name=Resolved_Name_L[Obs_ID_Inx]
        return Resolved_Name
def GC_Query(ObsID):
    #print "ObsID: ", ObsID
    Gname=Gname_Query(ObsID)
    #print "Gname: ", Gname
    #print "ObsID: ", ObsID
    if(Gname=="Error"):
        print "Error Gname: ", Gname
        #return "Error","Error"
        return np.nan,np.nan
    try:
        #print "NED Gname: ", Gname
        G_Data= Ned.query_object(Gname) #G_Data:-astropy.table.table.Table, Galaxy_Data, The Galaxy Data Table queried from NED
    except:
        Gname=Gname_Query(ObsID,NED_Resolved_Bool=True)
        #print "NED Gname: ", Gname
        G_Data= Ned.query_object(Gname) #G_Data:-astropy.table.table.Table, Galaxy_Data, The Galaxy Data Table queried from NED
    #print G_Data
    try:
        raGC=float(G_Data['RA(deg)'])
        decGC=float(G_Data['DEC(deg)'])
    except:
        raGC=float(G_Data['RA'])
        decGC=float(G_Data['DEC'])
    return raGC, decGC

def D25_Query(ObsID):
    ##Gname=Gname_Query(ObsID)
    #Evt2_File_H_L=File_Query_Code_5.File_Query(Gname,"evt2")
    #print "Evt2_File_H_L : ", Evt2_File_H_L
    #Evt2_File_L=Evt2_File_H_L[0]
    #Obs_ID=Evt2_File_L[0]
    Obs_ID=int(ObsID)
    #Evt2filepath=Evt2_File_Query(ObsID)
    #print "Obs_ID : ", Obs_ID
    #path_Source_Flux_Table=os.path.realpath('../SQL_Standard_File/Source_Flux_Table.csv') #Reltive Path
    path_Source_Flux_Table=os.path.realpath('/Volumes/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_Table.csv') #Absolute Path
    #print "path_Source_Flux_Table : ", path_Source_Flux_Table
    Source_Flux_Table=ascii.read(path_Source_Flux_Table)
    #print "Source_Flux_Table : \n", Source_Flux_Table
    Obs_ID_A=Source_Flux_Table["OBSID"]
    Obs_ID_L=list(Obs_ID_A)
    #print "Obs_ID_L : ", Obs_ID_L
    Obs_ID_Inx=Obs_ID_L.index(Obs_ID)
    #print "Obs_ID_Inx : ", Obs_ID_Inx
    Resolved_Name_A=Source_Flux_Table["resolvedObject"]
    #print "Resolved_Name_A : ", Resolved_Name_A
    Resolved_Name_L=list(Resolved_Name_A)
    #print "Resolved_Name_L : ", Resolved_Name_L
    Resolved_Name=Resolved_Name_L[Obs_ID_Inx]
    #print "Resolved_Name : ", Resolved_Name
    #Dia_Table = Ned.get_table(Gname, table='diameters') #Dia_Table:-astropy.table.table.Table, Diameter_Table, The Data table queried from NED that contains the infomation about the Major Axis of the input Galaxy Name
    Dia_Table = Ned.get_table(Resolved_Name, table='diameters') #Dia_Table:-astropy.table.table.Table, Diameter_Table, The Data table queried from NED that contains the infomation about the Major Axis of the input Galaxy Name
    #print type(Dia_Table)
    #print G_Data
    #print Dia_Table
    #print Dia_Table.colnames
    #print Dia_Table.meta
    #print Dia_Table.columns
    Dia_Table_Feq=Dia_Table['Frequency targeted'] #Dia_Table_Feq:-astropy.table.column.MaskedColumn, Diameter_Table_Fequency, The Array containing all named frequencies of light that are being used for the Major Axis Measurement
    #print Dia_Table['NED Frequency']
    #print Dia_Table_Feq
    #print type(Dia_Table_Feq)
    Dia_Table_Feq_L=list(Dia_Table_Feq) #Dia_Table_Feq_L:-List, Diameter_Table_Fequency_List, The list containing all named frequencies of light that are being used for the Major Axis Measurement
    #print Dia_Table_Feq_L
    Dia_Table_Num=Dia_Table['No.'] #Dia_Table_Num:-astropy.table.column.MaskedColumn, Diameter_Table_Number, The number Ned assigns to
    #print Dia_Table_Num
    #print type(Dia_Table_Num)
    Dia_Table_Num_L=list(Dia_Table_Num)
    #print Dia_Table_Num_L
    for i in range(0,len(Dia_Table_Feq_L)-1): #There is a bug here with index matching, The matched index isn't that same index for the major axis
        Cur_Feq=Dia_Table_Feq_L[i]
        #print Cur_Feq
        if(Cur_Feq=="RC3 D_25, R_25 (blue)"):
            Match_inx=i
            Match_Feq=Dia_Table_Feq_L[Match_inx]
            Match_Num=Dia_Table_Num_L[Match_inx]
            #Match_Num
            #print "Match_Feq ", Match_Feq
            #print "Match_inx ", Match_inx
            #print "Match_Num ", Match_Num
    #Dia_Table_Maj=Dia_Table['Major Axis']
    Dia_Table_Maj=Dia_Table['NED Major Axis']
    #print Dia_Table_Maj
    Dia_Table_Maj_L=list(Dia_Table_Maj)
    #print Dia_Table_Maj_L
    Dia_Table_Maj_Units=Dia_Table['Major Axis Unit']
    #print Dia_Table_Maj_Units
    Dia_Table_Maj_Units_L=list(Dia_Table_Maj_Units)
    #print Dia_Table_Maj_Units_L
    #print "i ", i
    D25_Maj=Dia_Table_Maj_L[Match_inx]
    #print "D25_Maj ", D25_Maj
    D25_Units=Dia_Table_Maj_Units[Match_inx]
    #print "D25_Units ", D25_Units
    #print type(Dia_Table)
    #print Dia_Table.info()
    #Dia_Table_2=Dia_Table[6]
    #print Dia_Table_2
    #Maj=Dia_Table_2[18]
    #print "Maj, ! ! !", Maj
    D25_S_Maj=D25_Maj/2.0
    D25_S_Maj_Deg=D25_S_Maj/3600.0
    D25_S_Maj_Arcmin=D25_S_Maj_Deg*60.0
    #return D25_S_Maj_Deg
    return D25_S_Maj_Arcmin

def Exposure_Time_Calc(ObsID):
    """
    #/Volumes/xray/simon/all_chandra_observations/12095/primary/
    query_path='/Volumes/xray/simon/all_chandra_observations/'+str(ObsID)+'/primary/*evt2*'
    evtfpath_L=glob.glob(query_path)
    if(len(evtfpath_L)>1):
        print str(ObsID)+" has "+str(len(evtfpath_L)-1)+"Additional Evt2 Files ! ! !"
        return "Error"
    evtfpath=evtfpath_L[0]
    """
    evtfpath=Evt2_File_Query(ObsID)
    #print "evtfpath: ", evtfpath
    #print "type(evtfpath): ", type(evtfpath)
    hdulist = fits.open(evtfpath)
    Exposure_Time=hdulist[1].header['EXPOSURE']
    #print "Exposure_Time: ", Exposure_Time
    #print "type(Exposure_Time): ", type(Exposure_Time)
    #/Volumes/xray/simon/all_chandra_observations/12095/primary/
    return Exposure_Time
def Source_Known_Flux_Finder(ObsID,Source_Num):
    Nearest_Neighbor_Hybrid_Coords_Fpath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    Chip_ID=Hybrid_Sources_Data.iloc[Source_Num-1]["Chip_ID"] #This should be Chip_ID
    #print Chip_ID
    evtfpath=Evt2_File_Query(ObsID)
    hdulist = fits.open(evtfpath)
    Obs_Date_Str=hdulist[1].header['DATE-OBS']
    #print Obs_Date_Str
    #print type(Obs_Date_Str)
    Obs_Date_L=Obs_Date_Str.split("-")
    #print Obs_Date_L
    Obs_Year_Str=Obs_Date_L[0]
    Obs_Year=int(Obs_Year_Str)
    #PIMMS_A=pd.read_csv('~/Desktop/PIMMS/'+PIMMSfname) #This Works but the path should be changed to the one in Research_Git, Not the Desktop Path
    #PIMMS_filepath=os.path.realpath('../PIMMS/'+PIMMSfname)
    PIMMS_filepath="/Volumes/xray/anthony/Research_Git/PIMMS/PIMMS_Data.csv"
    #print "PIMMS_filepath : ",PIMMS_filepath
    PIMMS_A=pd.read_csv(PIMMS_filepath)
    #print PIMMS_A
    Cycle_A=PIMMS_A['Cycle']
    #print Cycle_A
    Year_A=PIMMS_A['Year']
    #print Year_A
    ACIS_I_A=PIMMS_A['ACIS-I(erg/cm**2/s)']
    #print ACIS_I_A
    ACIS_S_A=PIMMS_A['ACIS-S(erg/cm**2/s)']
    #print ACIS_S_A
    Year_L=list(Year_A)
    #print Year_L
    #Obs_Year=2006 #This is a test value
    for i in range(0,len(Year_L)):
        #print i
        if(Year_L[i]==Obs_Year):
            Row_Inx=i
            break
    #print Row_Inx
    #if((Chip_ID==3) or (Chip_ID==7)): #For Front illuminated?
    #Chip_ID=0 #This is a test value
    if(Chip_ID<4):
        Flux_A=ACIS_I_A
    if(Chip_ID>=4):
        Flux_A=ACIS_S_A
    #print Flux_A
    Flux=Flux_A[Row_Inx]
    #print "Flux: ", Flux
    return Flux

def Limiting_Flux_Calc(ObsID,Source_Num):
    Nearest_Neighbor_Hybrid_Coords_Fpath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    Offaxis_Angle=Hybrid_Sources_Data.iloc[Source_Num-1]["Offaxis_Angle"]
    Offaxis_Angle_Floor=np.floor(Offaxis_Angle)
    #/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/NGC_3631/Flux_90_Files/3951/NGC_3631_ObsID_3951_Flux_90.txt
    Flux_90_Fpath_L=glob.glob("/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/*/Flux_90_Files/"+str(ObsID)+"/*_Flux_90.txt")
    if(len(Flux_90_Fpath_L)==0):
        return
    #try:
    Flux_90_Fpath=Flux_90_Fpath_L[0]
    #except:
        #print "ObsID: "+str(ObsID)+" Source: "+str(Source_Num)+" has an error"
        #print "Flux_90_Fpath_L: ", Flux_90_Fpath_L
    Flux_90_File=open(Flux_90_Fpath)
    Flux_90_Str=Flux_90_File.read()
    if("Background not in range of 4D Graph Data interpolation" in Flux_90_Str):
        #return 1.0
        return
    Flux_90_Str_L=Flux_90_Str.split("|")[1].split("\n")
    #print "Flux_90_Str_L Before: ", Flux_90_Str_L
    Flux_90_Str_L.pop(0)
    Flux_90_Str_L.pop(len(Flux_90_Str_L)-1)
    #print "Flux_90_Str_L: ", Flux_90_Str_L
    Flux_90_L=[]
    for Flux_90_Str in Flux_90_Str_L:
        Flux_90=float(Flux_90_Str)
        Flux_90_L.append(Flux_90)
    Flux_90_File.close()
    #print "Flux_90_L: ", Flux_90_L
    #print "len(Flux_90_L): ", len(Flux_90_L)
    if((Offaxis_Angle_Floor>-1) and (Offaxis_Angle_Floor<10)):
        Offaxis_Angle_Floor_Int=int(Offaxis_Angle_Floor)
        #try:
        Limiting_Flux=Flux_90_L[Offaxis_Angle_Floor_Int]
        #except:
            #print "ObsID: "+str(ObsID)+" Source: "+str(Source_Num)+" has an error"
            #print "Offaxis_Angle: ", Offaxis_Angle
            #print "Offaxis_Angle_Floor_Int: ", Offaxis_Angle_Floor_Int
        #print "type(Limiting_Flux): ",type(Limiting_Flux)
        #Limiting_Flux=np.fromstring(Limiting_Flux) #numpy.float64
        Limiting_Flux=np.float64(Limiting_Flux)
        return Limiting_Flux
    else:
        ##return "Outside_FOV" #For when the soruces is outside the reasonable field of veiw
        #return 1.0
        return
def Flux_Cut_Bool_Calc(Flux,Limiting_Flux):
    if((Limiting_Flux=="Outside_FOV") or (Limiting_Flux is None)):
        ##return "Outside_FOV"
        #return False
        #return 1.0
        return
    else:
        Flux_Cut_Bool=(Flux>Limiting_Flux) #The source is kept in the sample if True and removed from the sample if false
        return Flux_Cut_Bool
def Inside_FOV_Bool_Calc(ObsID,Source_Num):
    Nearest_Neighbor_Hybrid_Coords_Fpath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    Offaxis_Angle=Hybrid_Sources_Data.iloc[Source_Num-1]["Offaxis_Angle"]
    Offaxis_Angle_Floor=np.floor(Offaxis_Angle)
    if((Offaxis_Angle_Floor>-1) and (Offaxis_Angle_Floor<10)):
        return True
    else:
        return False
def Offaxis_Angle_Calc(ObsID,Source_Num):
    Nearest_Neighbor_Hybrid_Coords_Fpath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    Offaxis_Angle=Hybrid_Sources_Data.iloc[Source_Num-1]["Offaxis_Angle"]
    return Offaxis_Angle
def RA_DEC_Calc(ObsID,Source_Num):
    Nearest_Neighbor_Hybrid_Coords_Fpath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    RA=Hybrid_Sources_Data.iloc[Source_Num-1]["RA"]
    Dec=Hybrid_Sources_Data.iloc[Source_Num-1]["DEC"]
    return RA, Dec
def Distance_From_GC_Calc(ObsID,Source_Num):
    Source_RA,Source_DEC=RA_DEC_Calc(ObsID,Source_Num)
    Source_RA_Arcmin=Source_RA*60.0
    Source_DEC_Arcmin=Source_DEC*60.0
    GC_RA,GC_DEC=GC_Query(ObsID)
    #if(GC_RA=="Error"):
    if(GC_RA==np.nan):
        return
    GC_RA_Arcmin=GC_RA*60.0
    GC_DEC_Arcmin=GC_DEC*60.0
    dist=np.sqrt(((GC_DEC_Arcmin-Source_DEC_Arcmin)**2.0)+((GC_RA_Arcmin-Source_RA_Arcmin)**2.0))
    return dist
def Outside_D25_Bool_Calc(ObsID,Source_Num):
    D25_S_Maj_Arcmin=D25_Query(ObsID)
    Dist_Arcmin=Distance_From_GC_Calc(ObsID,Source_Num)
    Outside_D25_Bool=(float(Dist_Arcmin)>float(D25_S_Maj_Arcmin))
    return Outside_D25_Bool
def Distance_Galatic_Center_to_Aimpoint_Calc(ObsID):
    #pass
    GC_RA,GC_DEC=GC_Query(ObsID)
    #if(GC_RA=="Error"):
    if(GC_RA==np.nan):
        return
    #GC_RA_Arcmin=GC_RA*60.0
    #GC_DEC_Arcmin=GC_DEC*60.0
    Cur_Evt2_Filepath=Evt2_File_Query(ObsID)
    hdulist = fits.open(Cur_Evt2_Filepath)
    #Exposure_Time=hdulist[1].header['EXPOSURE']
    Pointing_RA=hdulist[1].header['RA_PNT']
    Pointing_Dec=hdulist[1].header['DEC_PNT']
    Pointing_Diff_RA=Pointing_RA-GC_RA
    Pointing_Diff_Dec=Pointing_Dec-GC_DEC
    Dist=np.sqrt((Pointing_Diff_RA**2.0)+(Pointing_Diff_Dec**2.0))
    Dist_Arcmin=Dist*60.0
    #Dist_L=[Cur_Evt2_ObsID,Dist_Arcmin]
    #Dist_H_L.append(Dist_L)
    #return Dist_H_L
    return Dist_Arcmin
def Source_Counts_To_Flux_Converter(fpath,Outfpath,CR_K=0.01):
    data=pd.read_csv(fpath)
    data['RAW_COUNTS[.3-7.5]']=data['RAW_COUNTS[.3-1]']+data['RAW_COUNTS[1-2.1]']+data['RAW_COUNTS[2.1-7.5]']
    data['BKG_COUNTS[.3-7.5]']=data['BKG_COUNTS[.3-1]']+data['BKG_COUNTS[1-2.1]']+data['BKG_COUNTS[2.1-7.5]']
    data['NET_COUNTS[.3-7.5]']=data['RAW_COUNTS[.3-7.5]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[.3-7.5]'])
    data['NET_COUNTS[.3-1]']=data['RAW_COUNTS[.3-1]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[.3-1]'])
    data['NET_COUNTS[1-2.1]']=data['RAW_COUNTS[1-2.1]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[1-2.1]'])
    data['NET_COUNTS[2.1-7.5]']=data['RAW_COUNTS[2.1-7.5]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[2.1-7.5]'])
    ObsID_A=data['OBSID']
    Exposure_Time_A=ObsID_A.apply(Exposure_Time_Calc) #This works but takes to long as it unnecessarily repeats getting the exposure time for every object rather than every ObsID
    #print "Exposure_Time_A:\n", Exposure_Time_A
    data['Exposure_Time']=Exposure_Time_A
    #data['Count_Rate[.3-7.5]']=data['RAW_COUNTS[.3-7.5]']/data['Exposure_Time'] #This should use net counts instead of raw counts!
    data['Count_Rate[.3-7.5]']=data['NET_COUNTS[.3-7.5]']/data['Exposure_Time']
    Src_A=data['SOURCE']
    #df['new_column'] = np.vectorize(fx)(df['A'], df['B'])
    data['Known_Flux[.3-7.5]']=np.vectorize(Source_Known_Flux_Finder)(ObsID_A,Src_A)
    #print "data:\n", data
    """
    print "data Grouped:\n", data.groupby("OBSID").groups
    Data_Grouped=data.groupby("OBSID")
    for name, group in Data_Grouped:
        print(name)
        print(group)
    """
    #F_U=F_K*((float(C)/Exposure_Time)/CR_K) # F_U:-float, Flux Unknown, This converts the current count value to the current flux value by assuming a linear relationship bewteen counts and flux, F_U is the calculated flux
    data['Flux[.3-7.5]']=data['Known_Flux[.3-7.5]']*(data['Count_Rate[.3-7.5]']/CR_K)
    data['Limiting_Flux[.3-7.5]']=np.vectorize(Limiting_Flux_Calc)(ObsID_A,Src_A)
    #data['Flux_Cut_Bool[.3-7.5]']=np.where(data['Flux[.3-7.5]']>data['Flux_Cut_Bool[.3-7.5]'])
    Flux_A=data['Flux[.3-7.5]']
    #print "type(Flux_A[3]): ", type(Flux_A[3])
    Limiting_Flux_A=data['Limiting_Flux[.3-7.5]']
    #print "type(Limiting_Flux_A[3]): ", type(Limiting_Flux_A[3])
    data['Flux_Cut_Bool[.3-7.5]']=np.vectorize(Flux_Cut_Bool_Calc)(Flux_A,Limiting_Flux_A)
    data["Inside_FOV_Bool"]=np.vectorize(Inside_FOV_Bool_Calc)(ObsID_A,Src_A)
    Coords=np.vectorize(RA_DEC_Calc)(ObsID_A,Src_A)
    data["RA"]=Coords[0]
    data["DEC"]=Coords[1]
    #Offaxis_Angle_Calc(ObsID,Source_Num)
    data["Offaxis_Angle"]=np.vectorize(Offaxis_Angle_Calc)(ObsID_A,Src_A)
    data["Gname_Modifed"]=np.vectorize(Gname_Query)(ObsID_A)
    #GC_Query(ObsID)
    GC_Coords=np.vectorize(GC_Query)(ObsID_A)
    data["RA_GC"]=GC_Coords[0]
    data["DEC_GC"]=GC_Coords[1]
    #D25_Query(ObsID)
    data["D25"]=np.vectorize(D25_Query)(ObsID_A)
    #Distance_From_GC_Calc(ObsID,Source_Num)
    data["Source_Distance_From_GC"]=np.vectorize(Distance_From_GC_Calc)(ObsID_A,Src_A)
    #Outside_D25_Bool_Calc(ObsID,Source_Num)
    data["Outside_D25_Bool"]=np.vectorize(Outside_D25_Bool_Calc)(ObsID_A,Src_A)
    #Distance_Galatic_Center_to_Aimpoint_Calc(ObsID)
    data["Distance_GC_to_Aimpoint"]=np.vectorize(Distance_Galatic_Center_to_Aimpoint_Calc)(ObsID_A)
    #with pd.option_context('display.max_rows', None):  # more options can be specified also
        #print "data:\n", data
    print "data:\n", data
    #Outfname=fpath.split("/")[0].split(".")[0]+"_Flux_Calc.csv"
    #print "Outfname: ", Outfname
    #Outpath=fpath.split("/")[:-1]

    data.to_csv(Outfpath)
    #print "data.dtypes:\n", data.dtypes
    #print "data:\n", data.where(data['Flux_Cut_Bool[.3-7.5]']==True,True,False)
    #print "data:\n", data.loc(data['Flux_Cut_Bool[.3-7.5]']==True)
    #print "data:\n", data.loc(True)
    #print data['Flux_Cut_Bool[.3-7.5]'].loc(True)

def Limiting_Flux_Model(F_L,m,F_L_E):
    #F_L_Corrected=m*np.log10(F_L)+np.log10(F_L_E)
    ##F_L_Corrected=m*F_L+F_L_E
    #F_L_Corrected=m*(10.0**F_L)+(10.0**F_L_E)
    #F_L_Corrected=m*(F_L**F_L_E)
    ##F_L_Corrected=10.0**(m*np.log10(F_L)+np.log10(F_L_E))
    F_L_Corrected=10.0**((m*np.log10(F_L))+np.log10(F_L_E))
    #F_L_Corrected=10.0**(np.log10(F_L)+np.log10(F_L_E))

    return F_L_Corrected

#print Limiting_Flux_Model(5.18967806355e-16,0.22222222222222232,1.3666666666666667e-15)

def Limiting_Flux_Model_Flux_Cut_Calc(Data,m,F_L_E,F_L_Bounds):
    #Data_Slice=Data[(Data['Limiting_Flux[.3-7.5]']>F_L_Bounds[0]) and (Data['Limiting_Flux[.3-7.5]']<F_L_Bounds[1])]
    Data_Slice=Data[(Data['Limiting_Flux[.3-7.5]']>F_L_Bounds[0]) & (Data['Limiting_Flux[.3-7.5]']<F_L_Bounds[1])]
    Data_Slice=Data_Slice[Data_Slice['Flux[.3-7.5]']>0]
    #print "F_L_Bounds: ", F_L_Bounds
    #print "Data Slice:\n", Data_Slice['Limiting_Flux[.3-7.5]']
    Flux_A=Data_Slice['Flux[.3-7.5]']
    Limiting_Flux_A=Data_Slice['Limiting_Flux[.3-7.5]']
    #print "Data Slice ({},{}): {}".format(str(F_L_Bounds[0],str(F_L_Bounds[1]),str(Limiting_Flux_A))
    Data_Slice["Limiting_Flux_Model"]=np.vectorize(Limiting_Flux_Model)(Limiting_Flux_A,m,F_L_E)
    Limiting_Flux_Model_A=Data_Slice["Limiting_Flux_Model"]
    Data_Slice['Model_Flux_Cut_Bool[.3-7.5]']=np.vectorize(Flux_Cut_Bool_Calc)(Flux_A,Limiting_Flux_Model_A)
    Model_Flux_Cut_A=Data_Slice['Model_Flux_Cut_Bool[.3-7.5]']
    #print "Model_Flux_Cut_A:\n", Model_Flux_Cut_A
    ##Model_Flux_Cut_Fraction=np.mean(Model_Flux_Cut_A) #Note: I don't know if NaN values effect this
    Model_Flux_Cut_Fraction=np.nanmean(Model_Flux_Cut_A) #Note: This might fix the NaN value issue
    return Model_Flux_Cut_Fraction, Model_Flux_Cut_A

def Modified_Standard_Deviation(A,Mu=0.90):
    N=len(A)
    #Modified_SD=np.sqrt((np.sum((A-Mu)**2.0))/N)
    Modified_SD=np.sqrt((np.sum((A-Mu)**2.0))/N)
    return Modified_SD

def Limiting_Flux_Model_Fitting(Data,Slope_Bounds=[-2.0,2.0],Flux_Error_Bounds=[1e-16,1e-14],Steps=200):
    Slope_A=np.linspace(Slope_Bounds[0],Slope_Bounds[1],Steps)
    print "Slope_A: ",Slope_A
    ##Flux_Error_A=np.linspace(Flux_Error_Bounds[0],Flux_Error_Bounds[1],Steps)
    Flux_Error_A=np.geomspace(Flux_Error_Bounds[0],Flux_Error_Bounds[1],Steps)
    print "Flux_Error_A: ", Flux_Error_A
    F_L_Bounds_Full=[5e-17,2e-15]
    ##F_L_Bounds_A=np.linspace(F_L_Bounds_Full[0],F_L_Bounds_Full[1],10)
    ##F_L_Bounds_A=np.linspace(F_L_Bounds_Full[0],F_L_Bounds_Full[1],20)
    F_L_Bounds_A=np.geomspace(F_L_Bounds_Full[0],F_L_Bounds_Full[1],20)
    F_L_Bounds_HL=[]
    for i in range(len(F_L_Bounds_A)-1):
        Cur_Bounds=[F_L_Bounds_A[i],F_L_Bounds_A[i+1]]
        F_L_Bounds_HL.append(Cur_Bounds)
    print "F_L_Bounds_HL: ", F_L_Bounds_HL
    Min_Modified_Standard_Dev=100000000000000000
    for Slope in Slope_A:
        for Flux_Error in Flux_Error_A:
            Flux_Cut_Frac_L=[]
            Flux_Cut_Frac_HL=[]
            for F_L_Bounds in F_L_Bounds_HL:
                Cur_Slice_Flux_Cut=Limiting_Flux_Model_Flux_Cut_Calc(Data,Slope,Flux_Error,F_L_Bounds)
                Cur_Slice_Flux_Cut_Frac=Cur_Slice_Flux_Cut[0]
                Flux_Cut_Frac_L.append(Cur_Slice_Flux_Cut_Frac)
                Flux_Cut_Frac_HL.append(Cur_Slice_Flux_Cut)
            Flux_Cut_Frac_A=np.array(Flux_Cut_Frac_L)
            Flux_Cut_Frac_Modified_SD=Modified_Standard_Deviation(Flux_Cut_Frac_A)
            if(Flux_Cut_Frac_Modified_SD<Min_Modified_Standard_Dev):
                Min_Modified_Standard_Dev=Flux_Cut_Frac_Modified_SD
                Cur_Best_Fit_Parameters=[Slope,Flux_Error]
                Cur_Best_Flux_Cut_Frac_L=Flux_Cut_Frac_L
                Cur_Best_Flux_Cut_HL=Flux_Cut_Frac_HL
    #print "Cur_Best_Flux_Cut_HL:\n", Cur_Best_Flux_Cut_HL
    print "Cur_Best_Flux_Cut_Frac_L: ", Cur_Best_Flux_Cut_Frac_L
    return [Cur_Best_Fit_Parameters,Min_Modified_Standard_Dev]


def Data_Analysis(fpath,Outside_D25_Bool=False,Flux_Model_B=False,Output_File_Ext="pdf"):
    data=pd.read_csv(fpath)
    if(Outside_D25_Bool==True):
        #data=data.drop(data[data.Outside_D25_Bool==False].index, inplace=True)
        #new_dataframe = a_dataframe[a_dataframe.B <= 3]
        data=data[data.Outside_D25_Bool]
        data=data.reset_index(drop=True)
    print "data.columns:\n", data.columns
    #print "data:\n", data
    with pd.option_context('display.max_columns', None):  # more options can be specified also
        print "data:\n", data
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data)
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Slope_Bounds=[0,100.0],Flux_Error_Bounds=[0,5])
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Flux_Error_Bounds=[1e-16,20e-16])
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Flux_Error_Bounds=[1e-16,40e-16])
    if(Flux_Model_B):
        Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Flux_Error_Bounds=[1e-16,1e-13])
        print "Flux_Cut_Fit: ", Flux_Cut_Fit
    Test_Bool=False
    Flux_A=data['Flux[.3-7.5]']
    Flux_Avg=Flux_A.mean()
    print "Flux_Avg: ", Flux_Avg
    Flux_L=list(Flux_A)
    Limiting_Flux_A=data['Limiting_Flux[.3-7.5]']
    Limiting_Flux_Avg=Limiting_Flux_A.mean()
    print "Limiting_Flux_Avg: ", Limiting_Flux_Avg
    Limiting_Flux_L=list(Limiting_Flux_A)
    Flux_Cut_Bool_A=data['Flux_Cut_Bool[.3-7.5]']
    Flux_Cut_Bool_L=list(Flux_Cut_Bool_A)
    #for Flux_Cut_Bool in Flux_Cut_Bool_L:
    n=0
    for i in range(0,len(Flux_Cut_Bool_L)):
        Flux_Cut_Bool=Flux_Cut_Bool_L[i]
        if(Flux_Cut_Bool==False):
            #print "True"
            ##print i
            ##print "Flux: ",Flux_L[i]
            ##print "Limiting_Flux: ", Limiting_Flux_L[i]
            Test_Bool=True
            n=n+1
    print Test_Bool
    print n
    if(Flux_Model_B):
        Model_Limiting_Flux_A=Limiting_Flux_Model(Limiting_Flux_A,Flux_Cut_Fit[0][0],Flux_Cut_Fit[0][1])
        print "Model_Limiting_Flux_A: ", Model_Limiting_Flux_A
    #print "Flux_A[2351]: ", Flux_A[2351]
    #print "Limiting_Flux_A[2351]: ", Limiting_Flux_A[2351]
    #print "Model_Limiting_Flux_A[2351]: ", Model_Limiting_Flux_A[2351]
    #Outside_D25_Bool_A=data["Outside_D25_Bool"]
    #a_dataframe.drop(a_dataframe[a_dataframe.B > 3].index, inplace=True)
    #data.drop(data[data.Outside_D25_Bool==False].index, inplace=True)
    #plt.plot(Limiting_Flux_A,Flux_A,".")
    #"""
    plt.loglog(Limiting_Flux_A,Flux_A,".",label="Source_Flux")
    plt.loglog(Limiting_Flux_A,Limiting_Flux_A,label="Limiting_Flux")
    if(Flux_Model_B):
        plt.loglog(Limiting_Flux_A,Model_Limiting_Flux_A,label="Model_Limiting_Flux")
    #"""
    """
    plt.plot(Limiting_Flux_A,Flux_A,".",label="Source_Flux")
    plt.plot(Limiting_Flux_A,Limiting_Flux_A,label="Limiting_Flux")
    plt.plot(Limiting_Flux_A,Model_Limiting_Flux_A,'.',label="Model_Limiting_Flux")
    plt.ylim(0.0,2e-18)
    """
    #plt.ylim(0.0,1e-14)
    plt.xlabel("Limiting_Flux (erg/cm**2/s absorbed flux)")
    plt.ylabel("Flux (erg/cm**2/s absorbed flux)")
    plt.legend()
    #plt.show()
    if(Outside_D25_Bool):
        plt.savefig("Source_Flux_Vs_Limiting_Flux_Outside_D25."+Output_File_Ext)
    else:
        plt.savefig("Source_Flux_Vs_Limiting_Flux."+Output_File_Ext)
    plt.clf()
    #"""
    #plt.hist(Limiting_Flux_A)
    ##plt.hist(Flux_A,bins=1000,log=True,range=(0,10.0**(-14.0)))
    ##plt.hist(Limiting_Flux_A,bins=1000,log=True,range=(0,10.0**(-14.0)))
    plt.hist(Flux_A,bins=100,log=False,range=(0,10.0**(-14.0)),label="Source_Flux")
    plt.hist(Limiting_Flux_A,bins=100,log=False,range=(0,10.0**(-14.0)),label="Limiting_Flux")
    plt.xlabel("Flux (erg/cm**2/s absorbed flux)")
    plt.ylabel("Number of Fluxes")
    plt.legend()
    #plt.show()
    if(Outside_D25_Bool):
        plt.savefig("Source_Flux_Vs_Limiting_Flux_Histogram_Outside_D25."+Output_File_Ext)
    else:
        plt.savefig("Source_Flux_Vs_Limiting_Flux_Histogram."+Output_File_Ext)
    plt.clf()
    #For Log(N)-Log(S)
    plt.hist(Limiting_Flux_A,bins=100,log=True,cumulative=-1,histtype='step',range=(0,2.5E-15))
    #plt.hist(Limiting_Flux_A,bins=100,log=True,cumulative=False,histtype='step',range=(0,2.5E-15))
    plt.xlabel("Limiting Flux (erg/cm**2/s absorbed flux)")
    plt.ylabel("N(>S)")
    #plt.legend()
    #plt.show()
    if(Outside_D25_Bool):
        plt.savefig("LogN_LogS_Outside_D25."+Output_File_Ext)
    else:
        plt.savefig("LogN_LogS."+Output_File_Ext)
    plt.clf()
#Exposure_Time_Calc(12095)
#Gname_Query(12095)
#print Gname_Query(6869)
#print GC_Query(12095)
#print GC_Query(6869)
#5905
#print GC_Query(5905)
#9552
#print GC_Query(9552)
#print D25_Query(12095)
#Source_Known_Flux_Finder(12095,1)
#print Limiting_Flux_Calc(12095,1)
#print Limiting_Flux_Calc(12095,5)
#print Limiting_Flux_Calc(2197,5)
#print Limiting_Flux_Calc(6869,1)
#print Distance_From_GC_Calc(12095,5)
#print Distance_From_GC_Calc(6869,5)
#print Outside_D25_Bool_Calc(12095,5)
#print Distance_From_GC_Calc(12095,20)
#print Outside_D25_Bool_Calc(12095,20)
#print Distance_Galatic_Center_to_Aimpoint_Calc(12095)
#print Distance_Galatic_Center_to_Aimpoint_Calc(6869)
#Source_Counts_To_Flux_Converter("/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small.csv","/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv",)
#Source_Counts_To_Flux_Converter("/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info.csv","/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv")
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv')
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv',Outside_D25_Bool=True)
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv')
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True)
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True,Output_File_Ext="png")
Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True)
