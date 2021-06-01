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
import re
import pyregion
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
def Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num):
    Nearest_Neighbor_Hybrid_Coords_Fpath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    Offaxis_Angle=Hybrid_Sources_Data.iloc[Source_Num-1]["Offaxis_Angle"]
    Offaxis_Angle_Floor=np.floor(Offaxis_Angle)
    Offaxis_Angle_Annulus_Number=int(Offaxis_Angle_Floor)
    return Offaxis_Angle_Annulus_Number
def Offaxis_Angle_Annulus_CCD_Incompleteness_Calc(ObsID,Source_Num):
    #evtfpath_L=glob.glob(query_path)
    #/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/NGC_3631/Area_Lists/3951/NGC_3631_acisf03951N004_evt2_Area_List.txt
    Offaxis_Angle_Annulus_Number=Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num)
    if(Offaxis_Angle_Annulus_Number>9):
        return np.nan
    Gname=Gname_Query(ObsID)
    #Offaxis_Angle_Annulus_Number=Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num)
    Area_List_Query_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+str(Gname)+"/Area_Lists/"+str(ObsID)+"/*evt2_Area_List.txt"
    Area_List_Path_L=glob.glob(Area_List_Query_Path)
    if(len(Area_List_Path_L)!=1):
        #print str(ObsID)+" has "+str(len(evtfpath_L)-1)+"Additional Evt2 Files ! ! !"
        print str(ObsID)+" has "+str(len(Area_List_Path_L))+" Area List Files ! ! !"
        ##return "Error"
        return np.nan
    Area_List_Path=Area_List_Path_L[0]
    Area_List_Path_Data=pd.read_csv(Area_List_Path,names=["Area_Ratios"])
    #print "Area_List_Path_Data\n", Area_List_Path_Data
    #return Area_List_Path_Data
    #Area_Ratio=Area_List_Path_Data[Offaxis_Angle_Annulus_Number]
    #print "Offaxis_Angle_Annulus_Number", Offaxis_Angle_Annulus_Number
    Area_Ratio=Area_List_Path_Data.iloc[Offaxis_Angle_Annulus_Number]
    #print type(Area_Ratio)
    return Area_Ratio
def Offaxis_Angle_Annulus_Area_Calc(ObsID,Source_Num):
    W=(1.0/60.0)
    Offaxis_Angle_Annulus_Number=Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num)
    if(Offaxis_Angle_Annulus_Number>9):
        return np.nan
    Annulus_Area=np.pi*(W**2.0)*((2.0*Offaxis_Angle_Annulus_Number)+1.0)
    return Annulus_Area
def Observed_Offaxis_Angle_Annulus_Area_Calc(ObsID,Source_Num):
    Annulus_Area=Offaxis_Angle_Annulus_Area_Calc(ObsID,Source_Num)
    Area_Ratio=Offaxis_Angle_Annulus_CCD_Incompleteness_Calc(ObsID,Source_Num)
    Observed_Offaxis_Angle_Annulus_Area=Annulus_Area*Area_Ratio
    return Observed_Offaxis_Angle_Annulus_Area
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
def Aimpoint_Coords_Calc(ObsID):
    Cur_Evt2_Filepath=Evt2_File_Query(ObsID)
    hdulist = fits.open(Cur_Evt2_Filepath)
    #Exposure_Time=hdulist[1].header['EXPOSURE']
    Pointing_RA=hdulist[1].header['RA_PNT']
    Pointing_Dec=hdulist[1].header['DEC_PNT']
    return Pointing_RA, Pointing_Dec
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
def Overlapping_ObsID_Calc(Data,Dist_Threshold=2.0):
    #pass
    #Aimpoint_Coords=np.vectorize(Aimpoint_Coords_Calc)(ObsID)
    RA_Aimpoint_A=Data["RA_Aimpoint"]
    RA_Aimpoint_L=list(RA_Aimpoint_A)
    DEC_Aimpoint_A=Data["DEC_Aimpoint"]
    DEC_Aimpoint_L=list(DEC_Aimpoint_A)
    ObsID_A=Data['OBSID']
    ObsID_L=list(ObsID_A)
    #Aimpoint_Coords=data['OBSID',"RA_Aimpoint","DEC_Aimpoint"]
    Overlapping_ObsID_HL=[]
    Overlapping_ObsID_Bool_L=[]
    for i in range(0,len(ObsID_L)):
        #j=i+1
        Overlapping_ObsID_Bool=False
        ObsID=ObsID_L[i]
        RA_Aimpoint=RA_Aimpoint_L[i]
        DEC_Aimpoint=DEC_Aimpoint_L[i]
        Overlapping_ObsID_L=[]
        #for j in range(i+1,len(ObsID_L)):
        for j in range(0,len(ObsID_L)):
            ObsID_Test=ObsID_L[j]
            if(ObsID_Test==ObsID):
                continue
            RA_Aimpoint_Test=RA_Aimpoint_L[j]
            DEC_Aimpoint_Test=DEC_Aimpoint_L[j]
            Aimpoint_Diff_RA=RA_Aimpoint-RA_Aimpoint_Test
            Aimpoint_Diff_DEC=DEC_Aimpoint-DEC_Aimpoint_Test
            Dist=np.sqrt((Aimpoint_Diff_RA**2.0)+(Aimpoint_Diff_DEC**2.0))
            if(Dist<Dist_Threshold):
                Overlapping_ObsID_Bool=True
                #print ObsID_Test
                Overlapping_ObsID_L.append(ObsID_Test)
        Overlapping_ObsID_Bool_L.append(Overlapping_ObsID_Bool)
        Overlapping_ObsID_L=list(set(Overlapping_ObsID_L))
        Overlapping_ObsID_L_Str=str(Overlapping_ObsID_L)
        Overlapping_ObsID_L_Str=Overlapping_ObsID_L_Str.replace(",",";")
        Overlapping_ObsID_HL.append(Overlapping_ObsID_L_Str)
    #Overlapping_ObsID_DF=pd.DataFrame(Overlapping_ObsID_HL)
    Data["Close_ObsIDs_Bool"]=Overlapping_ObsID_Bool_L
    Data["Close_ObsIDs"]=Overlapping_ObsID_HL
    return Data
def Duplicate_Source_Calc(Data,Dist_Threshold=(2.0/3600.0)): #Threshold needs to be determed. Temperary value used.
    #pass
    #Aimpoint_Coords=np.vectorize(Aimpoint_Coords_Calc)(ObsID)
    ObsID_A=Data['OBSID']
    ObsID_L=list(ObsID_A)
    Source_Num_A=Data['SOURCE']
    Source_Num_L=list(Source_Num_A)
    Source_RA_A=Data['RA']
    Source_RA_L=list(Source_RA_A)
    Source_Dec_A=Data['DEC']
    Source_Dec_L=list(Source_Dec_A)
    #Close_ObsIDs_Bool
    Close_ObsIDs_Bool_A=Data['Close_ObsIDs_Bool']
    Close_ObsIDs_Bool_L=list(Close_ObsIDs_Bool_A)
    #Close_ObsIDs
    Close_ObsIDs_A=Data['Close_ObsIDs']
    Close_ObsIDs_HL=list(Close_ObsIDs_A)
    #Aimpoint_Coords=data['OBSID',"RA_Aimpoint","DEC_Aimpoint"]
    Duplicate_Source_HL=[]
    Duplicate_Source_Bool_L=[]
    for i in range(0,len(ObsID_L)):
        Duplicate_Source_Bool=False
        Duplicate_Source_L=[]
        Close_ObsIDs_Bool=Close_ObsIDs_Bool_L[i]
        if(Close_ObsIDs_Bool==False):
            #Duplicate_Source_Bool=False
            Duplicate_Source_L=[]
            Duplicate_Source_Bool_L.append(Duplicate_Source_Bool)
            Duplicate_Source_HL.append(Duplicate_Source_L)
            continue
        ObsID=ObsID_L[i]
        Source_Num=Source_Num_L[i]
        Source_RA=Source_RA_L[i]
        Source_Dec=Source_Dec_L[i]
        #=re.split("[(),]",Cur_Raytrace_Reg)
        Close_ObsIDs_L_Str=Close_ObsIDs_HL[i]
        Close_ObsIDs_L_Strings=re.split("[\[\];]",Close_ObsIDs_L_Str) #Note: I am not sure if this is the correct regex
        print "Close_ObsIDs_L_Strings: ", Close_ObsIDs_L_Strings
        Close_ObsIDs_L=[]
        for ObsID_Str in Close_ObsIDs_L_Strings:
            if(ObsID_Str==""):
                continue
            ObsID_Int=int(ObsID_Str)
            Close_ObsIDs_L.append(ObsID_Int)
        #df.loc[df['column_name'].isin(some_values)]
        Data_Test=Data.loc[Data['OBSID'].isin(Close_ObsIDs_L)]
        ObsID_A_Test=Data_Test['OBSID']
        ObsID_L_Test=list(ObsID_A_Test)
        Source_Num_A_Test=Data_Test['SOURCE']
        Source_Num_L_Test=list(Source_Num_A_Test)
        Source_RA_A_Test=Data_Test['RA']
        Source_RA_L_Test=list(Source_RA_A_Test)
        Source_Dec_A_Test=Data['DEC']
        Source_Dec_L_Test=list(Source_Dec_A_Test)
        for j in range(0,len(ObsID_L_Test)):
            ObsID_Test=ObsID_L_Test[j]
            if(ObsID_Test==ObsID):
                continue
            Source_Num_Test=Source_Num_L_Test[j]
            Source_RA_Test=Source_RA_L_Test[j]
            Source_Dec_Test=Source_Dec_L_Test[j]
            Source_RA_Diff=Source_RA-Source_RA_Test
            Source_Dec_Diff=Source_Dec-Source_Dec_Test
            Dist=np.sqrt((Source_RA_Diff**2.0)+(Source_Dec_Diff**2.0))
            if(Dist<Dist_Threshold):
                Duplicate_Source_Bool=True
                #print ObsID_Test
                Duplicate_Source=[ObsID_Test,Source_Num_Test]
                Duplicate_Source_L.append(Duplicate_Source)
        Duplicate_Source_Bool_L.append(Duplicate_Source_Bool)
        Duplicate_Source_HL.append(Duplicate_Source_L)
        #Overlapping_ObsID_Bool_L.append(Overlapping_ObsID_Bool)
        #Overlapping_ObsID_L=list(set(Overlapping_ObsID_L))
        #Overlapping_ObsID_L_Str=str(Overlapping_ObsID_L)
        #Overlapping_ObsID_L_Str=Overlapping_ObsID_L_Str.replace(",",";")
        #Overlapping_ObsID_HL.append(Overlapping_ObsID_L_Str)
    #Overlapping_ObsID_DF=pd.DataFrame(Overlapping_ObsID_HL)
    Data["Duplicate_Source_Bool"]=Duplicate_Source_Bool_L
    Data["Duplicate_Sources"]=Duplicate_Source_HL
    return Data
def CCD_Region_Query(ObsID):
    #Region_File_Query_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+str(Gname)+"/Area_Lists/"+str(ObsID)+"/*"+str(ObsID)+"*_CCD_Regions_simple_region_modifed_Code.txt"
    Region_File_Query_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/*/Area_Lists/"+str(ObsID)+"/*"+str(ObsID)+"*_CCD_Regions_simple_region_modifed_Code.txt"
    Region_File_Path_L=glob.glob(Region_File_Query_Path)
    if(len(Region_File_Path_L)!=1):
        #print str(ObsID)+" has "+str(len(evtfpath_L)-1)+"Additional Evt2 Files ! ! !"
        print str(ObsID)+" has "+str(len(Region_File_Path_L))+" Simple Region Files ! ! !"
        return "Error"
    Region_File_Path=Region_File_Path_L[0]
    return Region_File_Path
def Annulus_Shape_String_Generator(X,Y,Annulus_Number,W):
    Inner_Radius=W*Annulus_Number
    Outer_Radius=W*(Annulus_Number+1)
    #shape1 ='annulus(' + str(X_Phys) +','+ str(Y_Phys)+','+ str(cur_inner_r)+','+ str(cur_r)+')' #shape1:-str, shape1, The shape string of the current area annulus.
    Annulus_Shape='annulus(' + str(X) +','+ str(Y)+','+ str(Inner_Radius)+','+ str(Outer_Radius)+')'
    return Annulus_Shape
def Region_Coords_Converter(ObsID,Points_L,To_RA_Dec_Bool=False): #Note: This currently converts from SKY coordinates to RA and Dec but I think it might have to go in reverse to be used with the CXCRegion pakage.
    evtfpath=Evt2_File_Query(ObsID)
    Point_Converted_L=[]
    for Point in Points_L:
        X=Points_L[0]
        Y=Points_L[1]
        if(To_RA_Dec_Bool):
            dmcoords(infile=str(evtfpath),x=str(X), y=str(Y), option='sky', verbose=0, celfmt='deg') # Runs the dmcoords CIAO tool, which converts coordinates like CHIP_ID to SKY, the tool is now being used to convert the physical coordinates in pixels to RA and Dec in decimal degrees.
            RA=dmcoords.ra
            Dec=dmcoords.dec
            Point_Converted=[RA,Dec]
        else:
            RA=X
            Dec=Y
            dmcoords(infile=str(evtfpath),ra=str(RA), dec=str(Dec), option='cel', verbose=0, celfmt='deg') # Runs the dmcoords CIAO tool, which converts coordinates like CHIP_ID to SKY, the tool is now being used to convert the RA and Dec in decimal degrees to SKY coodinates in pixels.
            X_Phys=dmcoords.x #X_Phys:-float, X_Physical, The sky plane X pixel coordinate in units of pixels
            Y_Phys=dmcoords.y #Y_Phys:-float, Y_Physical, The sky plane Y pixel coordinate in units of pixels
            Point_Converted=[X_Phys,Y_Phys]
        Point_Converted_L.append(Point_Converted)
    return Point_Converted_L
def Region_Lengths_Converter(Lengths_L):
    #2.03252032520325 the converstion factor is 2.03252032520325pix/arcsec or 0.4920+-0.0001 arcsec/pix
    Lengths_A=np.array(Lengths_L)
    Lengths_Coverted_A=0.4920*Lengths_A
    Lengths_Coverted_L=list(Lengths_Coverted_A)
    return Lengths_Coverted_L
def Region_Parse(Region_Str,Num_Points=1,Fpath_Bool=False):
    #Reg_Data=pyregion.open(fpath) #Reg_Data:-pyregion.ShapeList, Region_Data, The Region Data that is obtainted when the pyregion module converts the Simple Region File for it's use
    if(Fpath_Bool):
        fpath=Region_Str
        Reg_Data=pyregion.open(fpath)
    else:
        Reg_Data=pyregion.parse(Region_Str)
    Shape_List=[] #Shape_List:-List, Shape_List, The list of all Cur_Shape strings for an observation
    Shape_Data_L=[] #Shape_Data_L:-List, Shape_Data_List, The list of all data needed to create a simple CCD region for every CCD polygon, represented as a high list where each high element is a list containing the data for a perticular CCD polygon, ie. This list is in the form [[[Midpoint_X_1,Midpoint_Y_1],D_1,Angle_1],[[Midpoint_X_2,Midpoint_Y_2],D_2,Angle_2],...[[Midpoint_X_n,Midpoint_Y_n],D_1,Angle_n]], where 1,2..n are differnt CCD polygons
    for Cur_Reg_Data in Reg_Data: # Selects the current CCD shape polygon
        Point_L=[] #Point_L, Point_list, The list of points that make up the current polygon
        Dimensions_L=[]
        Format=Cur_Reg_Data.coord_format #Format:-str, Format, The current coordinate format that the Cur_Reg_Data is in, This code will always use "physical" coordinates
        #print type(Format)
        Coords_L=Cur_Reg_Data.coord_list #Coords_L:-List, Coordinates_List, This is the list of coordinates that makes up the points that make the the current polygon in the form [X1,Y1,X2,Y2...,Xn,Yn] for finte n
        #print type(Coords_L)
        Shape=Cur_Reg_Data.name #Shape:-str, Shape, The string name of the current type of shape, all CCD Region Files will consist of only polygon shapes
        #print type(Shape)
        #print "Format ", Format
        #print "Shape ", Shape
        #print "Coords_L ", Coords_L
        #for i in range(0,len(Coords_L)-1,2): # Selects each x value in the Coords_L, ie for the list [X1,Y1,X2,Y2...,Xn,Yn], it only selects X1,X2...Xn for finte n
        for i in range(0,(Num_Points*2)-1,2): # Selects each x value in the Coords_L, ie for the list [X1,Y1,X2,Y2...,Xn,Yn], it only selects X1,X2...Xn for finte n
            Cur_Point=[] #Cur_Point:-List, Current_Point, The current point in the Coords_L, Points are in the form [X,Y] for example [X1,Y1]
            #print i
            Cur_Point.append(Coords_L[i]) # Adds the current X to the Current Point
            Cur_Point.append(Coords_L[i+1]) # Adds the current Y to the Current Point after the X
            Point_L.append(Cur_Point) # Appends the current point to the Point_List
        for i in range(Num_Points*2,len(Coords_L)):
            Cur_Dimension=Coords_L[i]
            Dimensions_L.append(Cur_Dimension)
        #print "Point_L ", Point_L
        Cur_Shape_Data=[Shape,Point_L,Dimensions_L]
        Shape_Data_L.append(Cur_Shape_Data)
    return Shape_Data_L
def Region_Union(Region_L):
    if(len(Region_L)==0):
        return
    if(len(Region_L)==1):
        return Region_L[0]
    Unioned_Region=Region_L[0]
    for i in range(1,len(Region_L)):
        Region=Region_L[i]
        Unioned_Region=Unioned_Region+Region
    return Unioned_Region
def Region_Intersection(Region_L):
    if(len(Region_L)==0):
        return
    if(len(Region_L)==1):
        return Region_L[0]
    Intersected_Region=Region_L[0]
    for i in range(1,len(Region_L)):
        Region=Region_L[i]
        Intersected_Region=Intersected_Region*Region
    return Intersected_Region
def Limiting_Flux_Annulus_Intersection(ObsID,X_Aimpoint,Y_Aimpoint,Annulus_Number):
    #pass
    CCD_Region_Fpath=CCD_Region_Query(ObsID)
    Annulus_Shape_Str=Annulus_Shape_String_Generator(X_Aimpoint,Y_Aimpoint,Annulus_Number,W=(1.0/60.0))
    Annulus_Shape_Data=Region_Parse(Annulus_Shape_Str)[0]
    print "Annulus_Shape_Data: ", Annulus_Shape_Data
    Annulus_Shape_Point_L=Annulus_Shape_Data[1]
    print "Annulus_Shape_Point_L: ", Annulus_Shape_Point_L
    Annulus_Shape_Radii_L=Annulus_Shape_Data[2]
    print "Annulus_Shape_Radii_L: ", Annulus_Shape_Radii_L
    Annulus_Shape_Point_L_Converted=Region_Coords_Converter(ObsID,Annulus_Shape_Point_L)
    Annulus_Shape_Radii_L_Converted=Region_Lengths_Converter(Annulus_Shape_Radii_L)
    Annulus_Shape_Point_Converted=Annulus_Shape_Point_L_Converted[0]
    X_Converted=Annulus_Shape_Point_Converted[0]
    Y_Converted=Annulus_Shape_Point_Converted[1]
    Inner_Radius_Converted=Annulus_Shape_Radii_L_Converted[0]
    Outer_Radius_Converted=Annulus_Shape_Radii_L_Converted[1]
    Annulus_Shape_Str_Converted='annulus(' + str(X_Converted) +','+ str(Y_Converted)+','+ str(Inner_Radius_Converted)+','+ str(Outer_Radius_Converted)+')'
    Annulus_Region=CXCRegion(Annulus_Shape_Str_Converted) #Note: Annulus_Region should be converted to SKY coordinates
    CCD_Region=CXCRegion(CCD_Region_Fpath)
    Annulus_Intersection_Region=Region_Intersection([CCD_Region,Annulus_Region])
    return Annulus_Intersection_Region
def Source_Counts_To_Flux_Converter(fpath,Outfpath,CR_K=0.01):
    data=pd.read_csv(fpath)
    data['RAW_COUNTS[.3-7.5]']=data['RAW_COUNTS[.3-1]']+data['RAW_COUNTS[1-2.1]']+data['RAW_COUNTS[2.1-7.5]']
    data['BKG_COUNTS[.3-7.5]']=data['BKG_COUNTS[.3-1]']+data['BKG_COUNTS[1-2.1]']+data['BKG_COUNTS[2.1-7.5]']
    data['NET_COUNTS[.3-7.5]']=data['RAW_COUNTS[.3-7.5]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[.3-7.5]'])
    data['NET_COUNTS[.3-1]']=data['RAW_COUNTS[.3-1]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[.3-1]'])
    data['NET_COUNTS[1-2.1]']=data['RAW_COUNTS[1-2.1]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[1-2.1]'])
    data['NET_COUNTS[2.1-7.5]']=data['RAW_COUNTS[2.1-7.5]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[2.1-7.5]'])
    ObsID_A=data['OBSID']
    ##'''
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
    #Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num)
    data["Offaxis_Angle_Annulus_Number"]=np.vectorize(Offaxis_Angle_Annulus_Number_Calc)(ObsID_A,Src_A)
    data["Offaxis_Angle_Annulus_CCD_Incompleteness"]=np.vectorize(Offaxis_Angle_Annulus_CCD_Incompleteness_Calc)(ObsID_A,Src_A)
    data["Offaxis_Angle_Annulus_Area"]=np.vectorize(Offaxis_Angle_Annulus_Area_Calc)(ObsID_A,Src_A)
    data["Observed_Offaxis_Angle_Annulus_Area"]=np.vectorize(Observed_Offaxis_Angle_Annulus_Area_Calc)(ObsID_A,Src_A)
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
    ##'''
    #Aimpoint_Coords_Calc(ObsID)
    Aimpoint_Coords=np.vectorize(Aimpoint_Coords_Calc)(ObsID_A)
    data["RA_Aimpoint"]=Aimpoint_Coords[0]
    data["DEC_Aimpoint"]=Aimpoint_Coords[1]
    data["Distance_GC_to_Aimpoint"]=np.vectorize(Distance_Galatic_Center_to_Aimpoint_Calc)(ObsID_A)
    data=Overlapping_ObsID_Calc(data)
    data=Duplicate_Source_Calc(data)
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
def Limiting_Flux_Data_Slice_Calc(Data,F_L_Bounds):
    #Data_Slice=Data[(Data['Limiting_Flux[.3-7.5]']>F_L_Bounds[0]) and (Data['Limiting_Flux[.3-7.5]']<F_L_Bounds[1])]
    Data_Slice=Data[(Data['Limiting_Flux[.3-7.5]']>F_L_Bounds[0]) & (Data['Limiting_Flux[.3-7.5]']<F_L_Bounds[1])]
    Data_Slice=Data_Slice[Data_Slice['Flux[.3-7.5]']>0]
    #print "F_L_Bounds: ", F_L_Bounds
    #print "Data Slice:\n", Data_Slice['Limiting_Flux[.3-7.5]']
    return Data_Slice

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

def Limiting_Flux_Area_Calc(Data,F_L_Bounds):
    Data_Slice=Limiting_Flux_Data_Slice_Calc(Data,F_L_Bounds)
    Data_Slice=Data_Slice.drop_duplicates(subset=["OBSID","Offaxis_Angle_Annulus_Number"])
    Data_Slice=Data_Slice.reset_index(drop=True)
    print "F_L_Bounds: ", F_L_Bounds
    #with pd.option_context('display.max_columns', None):
    #print "Limiting_Flux_Area Data_Slice:\n", Data_Slice
    print "Limiting_Flux_Area Data_Slice:\n", Data_Slice[["OBSID","Gname_Modifed","RA","Offaxis_Angle_Annulus_Number"]]
    Observed_Offaxis_Angle_Annulus_Area_A=Data_Slice["Observed_Offaxis_Angle_Annulus_Area"]
    Observed_Offaxis_Angle_Annulus_Area_Sum=Observed_Offaxis_Angle_Annulus_Area_A.sum()
    return Observed_Offaxis_Angle_Annulus_Area_Sum

def Background_Source_Calc(Flux,Hardness_Str="S"):
    if(Hardness_Str=="S"):
        #Soft Equation
        N=370.0*((Flux/(2.0E-15))**(-0.85)) #In sources per deg^2
    if(Hardness_Str=="H"):
        #Hard Equation
        N=1200.0*((Flux/(2.0E-15))**(-1.0)) #In sources per deg^2
    return N

def Data_Analysis(fpath,Outside_D25_Bool=False,Flux_Model_B=False,Output_File_Ext="pdf"):
    data=pd.read_csv(fpath)
    if(Outside_D25_Bool==True):
        #data=data.drop(data[data.Outside_D25_Bool==False].index, inplace=True)
        #new_dataframe = a_dataframe[a_dataframe.B <= 3]
        data=data[data.Outside_D25_Bool]
        data=data.reset_index(drop=True)
    #Close_ObsIDs_Bool
    data=data[data.Close_ObsIDs_Bool==False]
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
    Limiting_Flux_Range=(7E-17,2.5E-15)
    #LogN_LogS_Hist_A=plt.hist(Limiting_Flux_A,bins=100,log=True,cumulative=-1,histtype='step',range=(0,2.5E-15))
    ##LogN_LogS_Hist_A=plt.hist(Limiting_Flux_A,bins=100,log=True,cumulative=-1,histtype='step',range=(0,2.5E-15))
    ##LogN_LogS_Hist_A=plt.hist(Limiting_Flux_A,bins=100,log=True,cumulative=-1,histtype='step',range=(7E-17,2.5E-15))
    LogN_LogS_Hist_A=plt.hist(Limiting_Flux_A,bins=100,log=True,cumulative=-1,histtype='step',range=Limiting_Flux_Range)
    #print "LogN_LogS_Hist_A:\n", LogN_LogS_Hist
    #LogN_LogS_Hist_Unsummed_A=plt.hist(Limiting_Flux_A,bins=100,log=True,histtype='step',range=(0,2.5E-15))
    LogN_LogS_Hist_Unsummed_A=plt.hist(Limiting_Flux_A,bins=100,log=True,histtype='step',range=Limiting_Flux_Range)
    plt.clf()
    LogN_LogS_Hist=LogN_LogS_Hist_A[0]
    print "LogN_LogS_Hist:\n", LogN_LogS_Hist
    LogN_LogS_Bins=LogN_LogS_Hist_A[1]
    print "LogN_LogS_Bins:\n", LogN_LogS_Bins
    print "type(LogN_LogS_Bins): ", type(LogN_LogS_Bins)
    Limiting_Flux_Area_L=[]
    for i in range(0,len(list(LogN_LogS_Bins))-1):
        j=i+1
        F_L_Bin_Bounds=[LogN_LogS_Bins[i],LogN_LogS_Bins[j]]
        Limiting_Flux_Area=Limiting_Flux_Area_Calc(data,F_L_Bin_Bounds)
        #print "Limiting_Flux_Area: ", Limiting_Flux_Area
        Limiting_Flux_Area_L.append(Limiting_Flux_Area)
    Limiting_Flux_Area_A=np.array(Limiting_Flux_Area_L)
    print "Limiting_Flux_Area_A:\n", Limiting_Flux_Area_A
    Limiting_Flux_Area_A_Summed=np.cumsum(Limiting_Flux_Area_A[::-1])[::-1]
    print "Limiting_Flux_Area_A_Summed:\n", Limiting_Flux_Area_A_Summed
    LogN_LogS_Source_Density_Hist_A=LogN_LogS_Hist/Limiting_Flux_Area_A_Summed
    LogN_LogS_Unsummed_Hist=LogN_LogS_Hist_Unsummed_A[0]
    LogN_LogS_Unsummed_Hist_Error_A=np.sqrt(LogN_LogS_Unsummed_Hist)
    #LogN_LogS_Sum=0
    LogN_LogS_Error_Sum=0.0
    LogN_LogS_Error_Sum_L=[]
    for LogN_LogS_Error in reversed(LogN_LogS_Unsummed_Hist_Error_A):
        LogN_LogS_Error_Sum=np.sqrt((LogN_LogS_Error_Sum**2.0)+(LogN_LogS_Error**2.0))
        LogN_LogS_Error_Sum_L.append(LogN_LogS_Error_Sum)
    LogN_LogS_Error_Sum_L=LogN_LogS_Error_Sum_L[::-1]
    #LogN_LogS_Source_Density_Error_A=np.sqrt(LogN_LogS_Hist)/Limiting_Flux_Area_A_Summed #Note: I'm not sure this is valid for errors since the array is summed
    LogN_LogS_Source_Density_Error_A=LogN_LogS_Error_Sum_L/Limiting_Flux_Area_A_Summed #Note: I'm not sure this is valid for errors since the array is summed
    print "LogN_LogS_Source_Density_Hist_A:\n", LogN_LogS_Source_Density_Hist_A
    LogN_LogS_Bins_A=np.array(LogN_LogS_Bins)
    LogN_LogS_Bins_Left_Edges_A=LogN_LogS_Bins_A[:-1]
    #plt.bar(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A)
    ##plt.plot(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A,label="Data")
    plt.semilogx(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A,label="Data",marker=".")
    plt.errorbar(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A,yerr=LogN_LogS_Source_Density_Error_A,linestyle="")
    ##plt.semilogx(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Hist/5.0,label="Test",marker=".")
    LogN_LogS_Soft_A=Background_Source_Calc(LogN_LogS_Bins_Left_Edges_A,Hardness_Str="S")
    #plt.plot(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Soft_A,label="Soft_Giacconi",linestyle="--")
    plt.semilogx(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Soft_A,label="Soft_Giacconi",linestyle="--")
    LogN_LogS_Hard_A=Background_Source_Calc(LogN_LogS_Bins_Left_Edges_A,Hardness_Str="H")
    print "LogN_LogS_Hard_A:\n", LogN_LogS_Hard_A
    #plt.plot(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Hard_A,label="Hard_Giacconi",linestyle="--")
    plt.semilogx(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Hard_A,label="Hard_Giacconi",linestyle="--")
    """
    Limiting_Flux_Paper_Soft_A=np.linspace(2.0E-16,3.0E-13,num=100)
    LogN_LogS_Soft_A=Background_Source_Calc(Limiting_Flux_Paper_Soft_A,Hardness_Str="S")
    plt.plot(Limiting_Flux_Paper_Soft_A,LogN_LogS_Soft_A)
    print "LogN_LogS_Soft_A:\n", LogN_LogS_Soft_A
    Limiting_Flux_Paper_Hard_A=np.linspace(2.0E-15,3.0E-14,num=100)
    LogN_LogS_Hard_A=Background_Source_Calc(Limiting_Flux_Paper_Hard_A,Hardness_Str="H")
    print "LogN_LogS_Hard_A:\n", LogN_LogS_Hard_A
    plt.plot(Limiting_Flux_Paper_Hard_A,LogN_LogS_Hard_A)
    """
    #plt.hist(Limiting_Flux_A,bins=100,log=True,cumulative=False,histtype='step',range=(0,2.5E-15))
    plt.ylim(0.0,2000)
    plt.xlim(0,2.5E-15)
    plt.xlabel("Limiting Flux (erg/cm**2/s absorbed flux)")
    plt.ylabel("N(>S)")
    plt.legend()
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
#print Offaxis_Angle_Annulus_CCD_Incompleteness_Calc(12095,1)
#print Offaxis_Angle_Annulus_CCD_Incompleteness_Calc(12095,5)
#print Offaxis_Angle_Annulus_Area_Calc(12095,5)
#print Offaxis_Angle_Annulus_Area_Calc(12095,1)
#print Observed_Offaxis_Angle_Annulus_Area_Calc(12095,5)
#print CCD_Region_Query("NGC_3631",3951)
#print CCD_Region_Query("NGC_3631",3)
#print CCD_Region_Query(3951)
#print Region_Parse("circle(0,10,100)")
#print Region_Parse("/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/NGC_3631/Area_Lists/3951/acisf03951_001N004_CCD_Regions_simple_region_modifed_Code.txt",Fpath_Bool=True)
#Limiting_Flux_Annulus_Intersection(0,0,5,CCD_Region="")
##Source_Counts_To_Flux_Converter("/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small.csv","/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv",)
#Source_Counts_To_Flux_Converter("/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info.csv","/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv")
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv')
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv',Outside_D25_Bool=True)
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv')
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True)
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True,Output_File_Ext="png")
##Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True)
