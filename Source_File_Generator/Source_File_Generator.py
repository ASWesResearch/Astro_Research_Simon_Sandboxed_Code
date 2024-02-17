import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
from astropy.io.fits import Header
from astropy.io import fits
import os
##from astroquery.ned import Ned
from astroquery.ipac.ned import Ned
import astropy.units as u
import sys
from astropy.io import ascii
import re
import pyregion
from ciao_contrib.runtool import *
from region import *
import urllib.request
import bs4 as bs
#def Evt2_File_Query(ObsID):
def File_Query(ObsID,ObsID_Path='/Volumes/expansion/ObsIDs/',key="evt2"):
    Query_Path=ObsID_Path+str(ObsID)+'/**/*'+str(key)+'*'
    fpath_L=glob.glob(Query_Path, recursive = True)
    print("fpath_L: ", fpath_L)
    if(len(fpath_L)==0):
        raise Exception("ObsID "+str(ObsID)+" has 0 "+str(key)+" Files ! ! !")
    if(len(fpath_L)!=1):
        ##raise Exception(str(ObsID)+" has "+str(len(fpath_L))+" "+str(key)+" Files ! ! !")
        for cur_fpath in fpath_L:
            if("repro" in cur_fpath):
                #print("Match")
                fpath=cur_fpath
                break
            elif("/new/" in cur_fpath):
                fpath=cur_fpath
                break
            else:
                raise Exception("ObsID "+str(ObsID)+" has "+str(len(fpath_L))+" "+str(key)+" Files ! ! !")
    else:
        fpath=fpath_L[0]
    return fpath
def Gname_Query(ObsID, Data_Path="/opt/xray/anthony/Research_Git/SQL_Standard_File/ocatResult_Modified.csv", Key="Gname"):
    #query_path='/Volumes/xray/simon/all_chandra_observations/'+str(ObsID)+'/primary/*evt2*'
    Glob_L=glob.glob(Data_Path)
    if(len(Glob_L)!=1):
        raise Exception("Data_Path has "+str(len(Glob_L))+" Matching Files ! ! !")
    Data_Path=Glob_L[0]
    Data=pd.read_csv(Data_Path)
    #print("Data: ", Data)
    ##Data_Reduced=Data[Data["Obs ID"].isin([ObsID])]
    Data_Reduced=Data[Data["Obs ID"] == ObsID]
    Data_Reduced=Data_Reduced.reset_index(drop=True)
    #print("Data_Reduced: ", Data_Reduced)
    Gname=Data_Reduced[Key].loc[0]
    #print("Gname: ", Gname)
    return Gname
def GC_Query(ObsID):
    #print "ObsID: ", ObsID
    Gname=Gname_Query(ObsID)
    #print("Gname: ", Gname)
    #print("type(Gname): ", type(Gname))
    #print("ObsID: ", ObsID)
    if(str(Gname)=="nan"):
        print("Error Gname: ", Gname)
        #return "Error","Error"
        return np.nan,np.nan
    try:
        #print "NED Gname: ", Gname
        G_Data= Ned.query_object(Gname) #G_Data:-astropy.table.table.Table, Galaxy_Data, The Galaxy Data Table queried from NED
    except:
        raise Exception("Galaxy name "+str(Gname)+" ObsID "+str(ObsID)+" Not Queryied from NED")
    #print G_Data
    try:
        raGC=float(G_Data['RA(deg)'])
        decGC=float(G_Data['DEC(deg)'])
    except:
        raGC=float(G_Data['RA'])
        decGC=float(G_Data['DEC'])
    return raGC, decGC
def Galaxy_Morph_Query(ObsID):
    #print "ObsID: ", ObsID
    Gname=Gname_Query(ObsID)
    if(isinstance(Gname, float)): #Checks for np.nan
        print("Error Gname: ", Gname)
        return np.nan
    #Gname_No_Whitespace=Gname.replace(" ","_")
    Gname_No_Whitespace=Gname.replace(" ","")
    #print("Gname_No_Whitespace: ", Gname_No_Whitespace)
    #print("Gname: ", Gname)
    #print("type(Gname): ", type(Gname))
    #print("ObsID: ", ObsID)
    if(str(Gname)=="nan"):
        print("Error Gname: ", Gname)
        #return "Error","Error"
        return np.nan
    #url="https://ned.ipac.caltech.edu/cgi-bin/NEDatt?objname=NGC_3631"
    url="https://ned.ipac.caltech.edu/cgi-bin/NEDatt?objname="+str(Gname_No_Whitespace)
    #print("url: ", url)
    source = urllib.request.urlopen(url).read()
    soup = bs.BeautifulSoup(source,'html.parser')
    soup.prettify()
    tables = soup.find_all("table")
    #print("len(tables): ", len(tables))
    if(len(tables)<2):
        print("Only "+str(len(tables))+" NED Morphology Tables Found!")
        print("tables\n:", tables)
        return np.nan
        #table=tables[0]
    #else:
    table=tables[1]
    #print("tables:\n", tables)
    #print("table:\n", table)
    #Row_L=[]
    table_rows = table.find_all('tr')
    Match_Found_Bool=False
    for tr in table_rows:
        td = tr.find_all('td')
        row = [i.text for i in td]
        #print(row)
        if(len(row)==8):
            #print("row[1]: ", row[1])
            if(row[1]=="1991RC3.9.C...0000d"):
                Match_Row=row
                Match_Found_Bool=True
        #Row_L.append(row)
    if(Match_Found_Bool==False):
        print(Gname+" Morphology not found")
        return np.nan
    Morph_Row=Match_Row
    #Morph_Row=Row_L[25]
    #print("Morph_Row: ", Morph_Row)
    Morph_Str=Morph_Row[0]
    return Morph_Str

def Morph_Reducer(ObsID):
    Morph_Str=Galaxy_Morph_Query(ObsID)
    if(isinstance(Morph_Str,float)):
        return np.nan
    Morph_L=["S","E","I"]
    for Morph in Morph_L:
        if(Morph in Morph_Str):
            return Morph

def D25_Query(ObsID):
    Gname=Gname_Query(ObsID)
    if(str(Gname)=="nan"):
        print("Error Gname: ", Gname)
        return np.nan,np.nan,np.nan
    print("ObsID: ", ObsID)
    print("Gname: ", Gname)
    try:
        Dia_Table = Ned.get_table(Gname, table='diameters') #Dia_Table:-astropy.table.table.Table, Diameter_Table, The Data table queried from NED that contains the infomation about the Major Axis of the input Galaxy Name
    except:
        return np.nan,np.nan,np.nan
    #Dia_Table = Ned.get_table(Resolved_Name, table='diameters') #Dia_Table:-astropy.table.table.Table, Diameter_Table, The Data table queried from NED that contains the infomation about the Major Axis of the input Galaxy Name
    #print(type(Dia_Table))
    #print(Dia_Table)
    #print(Dia_Table.colnames)
    #print(Dia_Table.meta)
    #print(Dia_Table.columns)
    Dia_Table_Feq=Dia_Table['Frequency targeted'] #Dia_Table_Feq:-astropy.table.column.MaskedColumn, Diameter_Table_Fequency, The Array containing all named frequencies of light that are being used for the Major Axis Measurement
    #print(Dia_Table['NED Frequency'])
    #print(Dia_Table_Feq)
    #print(type(Dia_Table_Feq))
    Dia_Table_Feq_L=list(Dia_Table_Feq) #Dia_Table_Feq_L:-List, Diameter_Table_Fequency_List, The list containing all named frequencies of light that are being used for the Major Axis Measurement
    #print(Dia_Table_Feq_L)
    Dia_Table_Num=Dia_Table['No.'] #Dia_Table_Num:-astropy.table.column.MaskedColumn, Diameter_Table_Number, The number Ned assigns to
    #print Dia_Table_Num
    #print type(Dia_Table_Num)
    Dia_Table_Num_L=list(Dia_Table_Num)
    #print Dia_Table_Num_L
    Match_Bool=False
    for i in range(0,len(Dia_Table_Feq_L)-1): #There is a bug here with index matching, The matched index isn't that same index for the major axis
        Cur_Feq=Dia_Table_Feq_L[i]
        #print Cur_Feq
        if(Cur_Feq=="RC3 D_25, R_25 (blue)"):
            Match_inx=i
            Match_Feq=Dia_Table_Feq_L[Match_inx]
            Match_Num=Dia_Table_Num_L[Match_inx]
            Match_Bool=True
            #Match_Num
            #print "Match_Feq ", Match_Feq
            #print "Match_inx ", Match_inx
            #print "Match_Num ", Match_Num
    if(Match_Bool==False):
        return np.nan,np.nan,np.nan
    #Dia_Table_Maj=Dia_Table['Major Axis']
    Dia_Table_Maj=Dia_Table['NED Major Axis']
    Dia_Table_Min=Dia_Table['NED Minor Axis']
    Dia_Table_Angle=Dia_Table['NED Position Angle']
    #print("Dia_Table_Angle:\n", Dia_Table_Angle)
    #print Dia_Table_Maj
    Dia_Table_Maj_L=list(Dia_Table_Maj)
    Dia_Table_Min_L=list(Dia_Table_Min)
    Dia_Table_Angle_L=list(Dia_Table_Angle)
    #print Dia_Table_Maj_L
    Dia_Table_Maj_Units=Dia_Table['Major Axis Unit']
    #print Dia_Table_Maj_Units
    Dia_Table_Maj_Units_L=list(Dia_Table_Maj_Units)
    #print Dia_Table_Maj_Units_L
    #print "i ", i
    D25_Maj=Dia_Table_Maj_L[Match_inx]
    D25_Min=Dia_Table_Min_L[Match_inx]
    D25_Angle=Dia_Table_Angle_L[Match_inx]
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
    D25_S_Min=D25_Min/2.0
    D25_S_Min_Deg=D25_S_Min/3600.0
    D25_S_Min_Arcmin=D25_S_Min_Deg*60.0
    #return D25_S_Maj_Deg
    ##return D25_S_Maj_Arcmin
    return (D25_S_Maj_Arcmin, D25_S_Min_Arcmin, D25_Angle)

def Galatic_Distance_Query(ObsID,Distance_Ref_Code="1988NBGC.C....0000T"):
    """
    Gname=Gname_Query(ObsID)
    if(str(Gname)=="nan"):
        print("Error Gname: ", Gname)
        return np.nan
    print("ObsID: ", ObsID)
    print("Gname: ", Gname)
    #result_table = Ned.query_object(Gname)
    #print("result_table:\n", result_table)
    #print("result_table keys: ", list(result_table.keys()))
    try:
        Dist_Table = Ned.get_table(Gname, table='redshifts') #Dia_Table:-astropy.table.table.Table, Diameter_Table, The Data table queried from NED that contains the infomation about the Major Axis of the input Galaxy Name
    except:
        return np.nan
    print("Dist_Table:\n", Dist_Table)
    print("Dist_Table keys: ", list(Dist_Table.keys()))
    """
    #print "ObsID: ", ObsID
    Gname=Gname_Query(ObsID)
    if(isinstance(Gname, float)): #Checks for np.nan
        print("Error Gname: ", Gname)
        return np.nan
    #Gname_No_Whitespace=Gname.replace(" ","_")
    Gname_No_Whitespace=Gname.replace(" ","")
    #print("Gname_No_Whitespace: ", Gname_No_Whitespace)
    #print("Gname: ", Gname)
    #print("type(Gname): ", type(Gname))
    #print("ObsID: ", ObsID)
    if(str(Gname)=="nan"):
        print("Error Gname: ", Gname)
        #return "Error","Error"
        return np.nan
    #url="https://ned.ipac.caltech.edu/cgi-bin/NEDatt?objname=NGC_3631"
    #url="https://ned.ipac.caltech.edu/cgi-bin/NEDatt?objname="+str(Gname_No_Whitespace)
    #url="http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=NGC1300"
    url="http://ned.ipac.caltech.edu/cgi-bin/nDistance?name="+str(Gname_No_Whitespace)
    print("url: ", url)
    source = urllib.request.urlopen(url).read()
    soup = bs.BeautifulSoup(source,'html.parser')
    soup.prettify()
    tables = soup.find_all("table")
    #print("len(tables): ", len(tables))
    if(len(tables)<2):
        print("Only "+str(len(tables))+" NED Distance Tables Found!")
        print("tables\n:", tables)
        return np.nan
        #table=tables[0]
    #else:
    table=tables[1]
    #print("tables:\n", tables)
    #print("table:\n", table)
    #Row_L=[]
    table_rows = table.find_all('tr')
    Match_Found_Bool=False
    for tr in table_rows:
        td = tr.find_all('td')
        row = [i.text for i in td]
        #print(row)
        #print(len(row))
        if(len(row)==10):
            #print("row[1]: ", row[1])
            #if(row[3]=="1997ApJS..109..333W"):
            #1988NBGC.C....0000T
            #if(row[4]=="1984A&AS...56..381B"):
            if(row[4]==Distance_Ref_Code):
                Match_Row=row
                Match_Found_Bool=True
        #Row_L.append(row)
    if(Match_Found_Bool==False):
        print(Gname+" Distance not found")
        return np.nan
    #Match_Row=table_rows[0] #For Testing
    Dist_Row=Match_Row
    #Dist_Row=Row_L[25]
    #print("Dist_Row: ", Dist_Row)
    Dist_Str=Dist_Row[2]
    Dist=float(Dist_Str)
    return Dist

def Luminosity_Calc(Flux,Distance):
    #Distance_cm=Distance*(3.086e+24)
    Distance_cm=Distance*(3.086*(10**24))
    #L=4.0*np.pi*(Distance**2.0)*Flux
    L=4.0*np.pi*(Distance_cm**2.0)*Flux
    return L

def Luminosity_Calc_Bulk(Data):
    Key_L=list(Data.keys())
    print("Key_L: ", Key_L)
    Flux_Key_L=[]
    for Key in Key_L:
        if "FLUX" in Key:
            Flux_Key_L.append(Key)
    print("Flux_Key_L: ", Flux_Key_L)
    for Flux_Key in Flux_Key_L:
        Flux_A=Data[Flux_Key]
        Dist_A=Data["Galatic_Distance"]
        Lum_Key=Flux_Key.replace("FLUX","LUM")
        #data[Lum_Key]=np.vectorize(Inside_FOV_Bool_Calc)(ObsID_A,Src_A)
        #Luminosity_Calc(Flux,Distance)
        Lum_A=np.vectorize(Luminosity_Calc)(Flux_A,Dist_A)
        Column_Index=int(Data.columns.get_loc("Exposure_Time"))
        #df.insert(loc = 0, column = 'Band', value = Band)
        Data.insert(loc = Column_Index, column = Lum_Key, value = Lum_A)
    return Data

def Start_Date_Calc(ObsID):
    Cur_Evt2_Filepath=File_Query(ObsID)
    hdulist = fits.open(Cur_Evt2_Filepath)
    Start_Date=hdulist[1].header['DATE-OBS']
    return Start_Date

def Exposure_Time_Calc(ObsID):
    evtfpath=File_Query(ObsID)
    hdulist = fits.open(evtfpath)
    Exposure_Time=hdulist[1].header['EXPOSURE']
    #print("Exposure_Time: ", Exposure_Time)
    return Exposure_Time

def Roll_Angle_Calc(ObsID):
    evtfpath=File_Query(ObsID)
    hdulist = fits.open(evtfpath)
    Roll_Angle=hdulist[1].header['ROLL_PNT']
    #print("Roll_Angle: ", Roll_Angle)
    return Roll_Angle

def Source_Known_Flux_Finder(ObsID,Source_Num,Effective_Area_Correction_Bool=True):
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
    PIMMS_A=pd.read_csv(PIMMS_filepath) #BUG: Potential Major Bug Here! If Simon's soure counts are already effective area corrected. Then using the PIMMS data again must do the correction twice erroneously leading to the wrong result!
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
    if(Effective_Area_Correction_Bool):
        for i in range(0,len(Year_L)):
            #print i
            if(Year_L[i]==Obs_Year):
                Row_Inx=i
                break
    if(Effective_Area_Correction_Bool==False):
        Row_Inx=3
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
def Inside_FOV_Bool_Calc(ObsID,Source_Num,Reg_Path="/opt/xray/anthony/expansion_backup/Hybrid_Regions/"):
    Nearest_Neighbor_Hybrid_Coords_Fpath=Reg_Path+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    Offaxis_Angle=Hybrid_Sources_Data.iloc[Source_Num-1]["Offaxis_Angle"]
    Offaxis_Angle_Floor=np.floor(Offaxis_Angle)
    if((Offaxis_Angle_Floor>-1) and (Offaxis_Angle_Floor<10)):
        return True
    else:
        return False
def Offaxis_Angle_Calc(ObsID,Source_Num,Reg_Path="/opt/xray/anthony/expansion_backup/Hybrid_Regions/"):
    Nearest_Neighbor_Hybrid_Coords_Fpath=Reg_Path+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    Offaxis_Angle=Hybrid_Sources_Data.iloc[Source_Num-1]["Offaxis_Angle"]
    return Offaxis_Angle
def Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num,Reg_Path="/opt/xray/anthony/expansion_backup/Hybrid_Regions/"):
    Nearest_Neighbor_Hybrid_Coords_Fpath=Reg_Path+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
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
    if(str(Gname)=="nan"):
        print("Error Gname: ", Gname)
        return np.nan
    #Offaxis_Angle_Annulus_Number=Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num)
    Area_List_Query_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+str(Gname)+"/Area_Lists/"+str(ObsID)+"/*evt2_Area_List.txt"
    Area_List_Path_L=glob.glob(Area_List_Query_Path)
    if(len(Area_List_Path_L)!=1):
        #print str(ObsID)+" has "+str(len(evtfpath_L)-1)+"Additional Evt2 Files ! ! !"
        print(str(ObsID)+" has "+str(len(Area_List_Path_L))+" Area List Files ! ! !")
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
def RA_DEC_Calc(ObsID,Source_Num,Reg_Path="/opt/xray/anthony/expansion_backup/Hybrid_Regions/"):
    Nearest_Neighbor_Hybrid_Coords_Fpath=Reg_Path+str(ObsID)+"/Nearest_Neighbor_Hybrid_Sources_ObsID_"+str(ObsID)+"_Coords.csv"
    Hybrid_Sources_Data=pd.read_csv(Nearest_Neighbor_Hybrid_Coords_Fpath)
    #print Hybrid_Sources_Data
    RA=Hybrid_Sources_Data.iloc[Source_Num-1]["RA"]
    Dec=Hybrid_Sources_Data.iloc[Source_Num-1]["DEC"]
    return RA, Dec

def Deg_to_Rad(Angle):
    Angle=(np.pi/180.0)*Angle
    return Angle

def Rad_to_Deg(Angle):
    Angle=(180.0/np.pi)*Angle
    return Angle

def Angle_Convert(Angle):
    """
    Angle:-float, an anlge in radians

    Converts angles to be in the range of 0 rad<New_Angle<2pi rad
    """
    Full_Cricle=2.0*np.pi
    #print("Full_Cricle: ", Full_Cricle)
    if(Angle<0):
        New_Angle=Angle+Full_Cricle
        #print("New_Angle After: ", New_Angle)
        return New_Angle
    if(Angle>Full_Cricle):
        New_Angle=Angle-Full_Cricle
        return New_Angle
    else:
    #if((Angle>0) and (Angle<Full_Cricle)):
        New_Angle=Angle
    print("New_Angle Out: ", New_Angle)
    return New_Angle

def Haversine_Distance(x1,x2,y1,y2):
    dx=x2-x1
    dy=y2-y1
    B=np.sqrt(((np.sin(dy/2.0))**2.0)+((np.cos(y1)*np.cos(y2))*((np.sin(dx/2.0))**2.0)))
    if(y1<0): #Note: This might not work for galaxies on the equator! This needs to be tested! #Note: This works on the equator as well!
        B=-1.0*B
    Have_Dist=2.0*np.arcsin(B)
    return Have_Dist

def Distance_From_GC_Calc(ObsID,Source_Num):
    Source_RA,Source_DEC=RA_DEC_Calc(ObsID,Source_Num)
    Source_RA_Arcmin=Source_RA*60.0
    Source_DEC_Arcmin=Source_DEC*60.0
    GC_RA,GC_DEC=GC_Query(ObsID)
    x1_Rad=Deg_to_Rad(GC_RA)
    x2_Rad=Deg_to_Rad(Source_RA)
    y1_Rad=Deg_to_Rad(GC_DEC)
    y2_Rad=Deg_to_Rad(Source_DEC)
    Have_Dist=Haversine_Distance(x1_Rad,x2_Rad,y1_Rad,y2_Rad)
    Have_Dist=np.abs(Have_Dist)
    Have_Dist_Deg=Rad_to_Deg(Have_Dist)
    Have_Dist_Arcmin=Have_Dist_Deg*60.0
    print("Have_Dist_Arcmin: ", Have_Dist_Arcmin)
    dx_Have=Haversine_Distance(x1_Rad,x2_Rad,y1_Rad,y1_Rad)
    #if(GC_RA=="Error"):
    if(GC_RA==np.nan):
        return
    GC_RA_Arcmin=GC_RA*60.0
    GC_DEC_Arcmin=GC_DEC*60.0
    dist=np.sqrt(((GC_DEC_Arcmin-Source_DEC_Arcmin)**2.0)+((GC_RA_Arcmin-Source_RA_Arcmin)**2.0))
    #return dist
    return Have_Dist_Arcmin

def Outside_D25_Bool_Calc(ObsID,Source_Num):
    #D25_S_Maj_Arcmin=D25_Query(ObsID)
    D25_Info=D25_Query(ObsID)
    D25_S_Maj_Arcmin=D25_Info[0]
    if(D25_S_Maj_Arcmin==np.nan):
        return np.nan
    Dist_Arcmin=Distance_From_GC_Calc(ObsID,Source_Num)
    Outside_D25_Bool=(float(Dist_Arcmin)>float(D25_S_Maj_Arcmin))
    return Outside_D25_Bool

def Phi_Query(ObsID, Source_Num, Data_Path="/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All.csv", Key="PHI"): #Note: It would be better to refernce another file
    #query_path='/Volumes/xray/simon/all_chandra_observations/'+str(ObsID)+'/primary/*evt2*'
    Glob_L=glob.glob(Data_Path)
    if(len(Glob_L)!=1):
        raise Exception("Data_Path has "+str(len(Glob_L))+" Matching Files ! ! !")
    Data_Path=Glob_L[0]
    Data=pd.read_csv(Data_Path)
    #print("Data: ", Data)
    ##Data_Reduced=Data[Data["Obs ID"].isin([ObsID])]
    Data_Reduced=Data[Data["ObsID"] == ObsID]
    Data_Reduced=Data_Reduced.reset_index(drop=True)
    #print("Data_Reduced: ", Data_Reduced)
    #Phi=Data_Reduced[Key].loc[0]
    Phi_A=Data_Reduced[Key]
    #print("Phi_A:\n", Phi_A)
    #Offaxis_Angle=Hybrid_Sources_Data.iloc[Source_Num-1]["Offaxis_Angle"]
    Phi=Phi_A.iloc[Source_Num-1]
    return Phi

def Psi_Calc(ObsID, Source_Num):
    Roll_Angle=Roll_Angle_Calc(ObsID)
    #print("Roll_Angle: ", Roll_Angle)
    Phi=Phi_Query(ObsID, Source_Num)
    #print("Phi: ", Phi)
    Psi=Roll_Angle-Phi+90
    return Psi

def Lambda_Calc(ObsID, Source_Num):
    Psi=Psi_Calc(ObsID, Source_Num)
    D25_Info=D25_Query(ObsID)
    Galaxy_Position_Angle=D25_Info[2]
    if(Galaxy_Position_Angle==np.nan):
        return np.nan
    #Lambda=Psi-Galaxy_Position_Angle-270
    Lambda=Psi+Galaxy_Position_Angle-180
    return Lambda

def Deg_to_Rad(Angle):
    Angle=(np.pi/180.0)*Angle
    return Angle

def Rad_to_Deg(Angle):
    Angle=(180.0/np.pi)*Angle
    return Angle

def Angle_Convert(Angle):
    """
    Angle:-float, an anlge in radians

    Converts angles to be in the range of 0 rad<New_Angle<2pi rad
    """
    Full_Cricle=2.0*np.pi
    #print("Full_Cricle: ", Full_Cricle)
    if(Angle<0):
        New_Angle=Angle+Full_Cricle
        #print("New_Angle After: ", New_Angle)
        return New_Angle
    if(Angle>Full_Cricle):
        New_Angle=Angle-Full_Cricle
        return New_Angle
    else:
    #if((Angle>0) and (Angle<Full_Cricle)):
        New_Angle=Angle
    print("New_Angle Out: ", New_Angle)
    return New_Angle

def Ellipse_Radius_Calc(a,b,Ang_Rel):
    """
    a:-float, Semi-major axis of source ellipse
    b:-float, Semi-minor axis of source ellipse
    Ang_Rel:-float, Polar angle at which the ellipse radius will be calculated in radians

    Calcuates the radius of an ellipse at a given polar angle
    """
    if((a==0) or (b==0)):
        return 0.0
    r=(a*b)/(np.sqrt(((b*np.cos(Ang_Rel))**2.0)+((a*np.sin(Ang_Rel))**2.0)))
    return r

def Source_Overlap_Calc(x1,x2,y1,y2,a1,a2,b1,b2,rot1,rot2,M=2.0):
    """
    x1:-float, X Postion of First source ellipse
    x2:-float, X Postion of Second source ellipse
    y1:-float, Y Postion of First source ellipse
    y2:-float, Y Postion of Second source ellipse
    a1:-float, Semi-major axis of First source ellipse
    a2:-float, Semi-major axis of Second source ellipse
    b1:-float, Semi-minor axis of First source ellipse
    b2:-float, Semi-minor axis of Second source ellipse
    rot1:-float, rotation angle of First source ellipse
    rot2:-float, rotation angle of Second source ellipse

    This function takes the coordinates for two source ellipse region shapes as an input and then returns a list of booleans where
    the first element bool is True when the first source has the second source in its backgournd area and the second element bool is True
    when the first soruce area is overlapping with the second source area.
    """
    dx=x2-x1
    print("x2: ", x2)
    print("x1: ", x1)
    #dx=7.730 #For Testing
    print("dx Before: ", dx)
    dy=y2-y1
    print("dy: ", dy)
    x1_Deg=x1/60.0
    x1_Rad=Deg_to_Rad(x1_Deg)
    x2_Deg=x2/60.0
    x2_Rad=Deg_to_Rad(x2_Deg)
    y1_Deg=y1/60.0
    y1_Rad=Deg_to_Rad(y1_Deg)
    y2_Deg=y2/60.0
    y2_Rad=Deg_to_Rad(y2_Deg)
    Have_Dist=Haversine_Distance(x1_Rad,x2_Rad,y1_Rad,y2_Rad)
    Have_Dist=np.abs(Have_Dist)
    Have_Dist_Deg=Rad_to_Deg(Have_Dist)
    Have_Dist_Arcmin=Have_Dist_Deg*60.0
    #print("Have_Dist: ", Have_Dist)
    dx_Have=Haversine_Distance(x1_Rad,x2_Rad,y1_Rad,y1_Rad)
    dx_Have_Deg=Rad_to_Deg(dx_Have)
    dx_Have_Arcmin=dx_Have_Deg*60.0
    dx=dx_Have_Arcmin
    print("dx After: ", dx)
    if((dx==0) and (dy==0)):
        #print "Same_Source"
        return [False,False]
    Position_Angle_1=np.arctan2(dy,dx)
    #print("Position_Angle_1 deg: ", Position_Angle_1*(180.0/np.pi))
    #print("rot1: ", rot1)
    Position_Angle_1=Angle_Convert(Position_Angle_1)
    print("Position_Angle_1 deg: ", Position_Angle_1*(180.0/np.pi))
    Position_Angle_2=Position_Angle_1+np.pi
    Position_Angle_2=Angle_Convert(Position_Angle_2)
    Dist=np.sqrt((dx**2.0)+(dy**2.0))
    #Dist=Dist*60.0 #For degrees to arcmin (not always used)
    rot1_Rad=Deg_to_Rad(rot1)
    rot1_Rad=Angle_Convert(rot1_Rad)
    rot2_Rad=Deg_to_Rad(rot2)
    rot2_Rad=Angle_Convert(rot2_Rad)
    #print("Position_Angle_1: ", Position_Angle_1)
    #print("rot1_Rad: ", rot1_Rad)
    print("rot1_Rad deg: ", rot1_Rad*(180.0/np.pi))
    Relative_Angle_1=Position_Angle_1-rot1_Rad #Need to make sure rot1 is in radians
    #print("Relative_Angle_1: ", Relative_Angle_1)
    Relative_Angle_2=Position_Angle_2-rot2_Rad #Need to make sure rot2 is in radians
    print("Relative_Angle_1 deg Before: ", Relative_Angle_1*(180.0/np.pi))
    Relative_Angle_1=Angle_Convert(Relative_Angle_1)
    print("Relative_Angle_1 deg After: ", Relative_Angle_1*(180.0/np.pi))
    Relative_Angle_2=Angle_Convert(Relative_Angle_2)
    r1=Ellipse_Radius_Calc(a1,b1,Relative_Angle_1)
    r1_Background=M*r1
    r2=Ellipse_Radius_Calc(a2,b2,Relative_Angle_2)
    r1_x=r1*np.cos(Relative_Angle_1)
    r1_x_Deg=r1_x/60.0
    r1_x_Rad=Deg_to_Rad(r1_x_Deg)
    r1_y=r1*np.sin(Relative_Angle_1)
    r1_y_Deg=r1_y/60.0
    r1_y_Rad=Deg_to_Rad(r1_y_Deg)
    r1_Have=Haversine_Distance(x2_Rad,r1_x_Rad,y1_Rad,r1_y_Rad)
    r1_Have=np.abs(r1_Have)
    Dgs=Dist-r1-r2
    Dgs_Have=Have_Dist-r1_Have-r2 #Note: r2 not converted yet.
    #print("Have_Dist_Arcmin: ", Have_Dist_Arcmin)
    #print("r1_Have: ", r1_Have)
    print("Dist: ", Dist)
    print("r1: ", r1)
    #print("r2: ", r2)
    print("Dgs: ", Dgs)
    #print("Dgs_Have: ", Dgs_Have)
    Dgb=Dist-r1_Background-r2
    #print "Dgb: ", Dgb
    Background_Overlap_Bool=(Dgb<0)
    Source_Overlap_Bool=(Dgs<0)
    Overlap_List=[Background_Overlap_Bool,Source_Overlap_Bool]
    return Overlap_List

def Outside_Elliptical_D25_Bool_Calc(ObsID,Source_Num):
    print("Source_Num: ", Source_Num)
    Source_RA,Source_DEC=RA_DEC_Calc(ObsID,Source_Num)
    print("Source Coords: ", "("+str(Source_RA)+","+str(Source_DEC)+")")
    Source_RA_Arcmin=Source_RA*60.0
    Source_DEC_Arcmin=Source_DEC*60.0
    print("Source Coords Arcmin: ", "("+str(Source_RA_Arcmin)+","+str(Source_DEC_Arcmin)+")")
    print("Source Coords Arcmin deg: ", "("+str(Source_RA_Arcmin/60.0)+","+str(Source_DEC_Arcmin/60.0)+")")
    GC_RA,GC_DEC=GC_Query(ObsID)
    print("GC Coords: ", "("+str(GC_RA)+","+str(GC_DEC)+")")
    #if(GC_RA=="Error"):
    if(GC_RA==np.nan):
        return
    GC_RA_Arcmin=GC_RA*60.0
    GC_DEC_Arcmin=GC_DEC*60.0
    print("GC Coords Arcmin: ", "("+str(GC_RA_Arcmin)+","+str(GC_DEC_Arcmin)+")")
    print("GC Coords Arcmin deg: ", "("+str(GC_RA_Arcmin/60.0)+","+str(GC_DEC_Arcmin/60.0)+")")

    D25_Info=D25_Query(ObsID)
    D25_Maj=D25_Info[0]
    if(D25_Maj==np.nan):
        return
    #print(isinstance(D25_Query(7087)[2],np.ma.core.MaskedConstant))
    D25_Min=D25_Info[1]
    D25_Angle=D25_Info[2] #Note: I think this needs to be converted to the correct angle measured from the +X-axis to the +Y-axis (W to N)
    print("D25_Angle: ", D25_Angle)
    if(isinstance(D25_Min,np.ma.core.MaskedConstant) or isinstance(D25_Angle,np.ma.core.MaskedConstant)): #In the event that the ellipse is assumed to be a cirlce by RC3
        print(str(ObsID)+" has a circular D25 region")
        return Outside_D25_Bool_Calc(ObsID,Source_Num)
    """
    #Beta=D25_Angle-90.0 #This is the converstion to the correct angle measured from the +X-axis to the +Y-axis (W to N)
    #Beta=(-1.0*D25_Angle)+90.0 #This is the converstion to the correct angle measured from the +X-axis to the +Y-axis (W to N)
    #Beta=(-1.0*D25_Angle)+90.0 #This is the converstion to the correct angle measured from the +X-axis to the +Y-axis (W to N)
    Beta=90.0-D25_Angle #This is the converstion to the correct angle measured from the +X-axis to the +Y-axis (W to N)
    #Beta=90.0+D25_Angle #This is the converstion to the correct angle measured from the +X-axis to the +Y-axis (W to N)
    #Beta=D25_Angle-270
    #if(Beta<0):
        #Beta=-1.0*Beta
        #Beta=Beta+360.0
    if(D25_Angle>90.0):
        Beta=D25_Angle-90.0
        #Beta=D25_Angle+90.0
        #Beta=Beta+180.0
    print("Beta: ", Beta)
    """
    #dx=Source_RA_Arcmin-GC_RA_Arcmin
    dx=GC_RA_Arcmin-Source_RA_Arcmin
    dy=Source_DEC_Arcmin-GC_DEC_Arcmin
    if((D25_Angle>90) and (dx>0) and (dy>0)):
        print("Q1 Condition")
        Theta=D25_Angle-90
    elif((D25_Angle<90) and (dx>0) and (dy<0)):
        print("Q4 Condition")
        Theta=D25_Angle+270
        #GC_RA_Arcmin=-1.0*GC_RA_Arcmin
        #Source_RA_Arcmin=-1.0*Source_RA_Arcmin
    else:
        Theta=D25_Angle+90
    print("Theta: ", Theta)

    #Source_Overlap_Calc(x1,x2,y1,y2,a1,a2,b1,b2,rot1,rot2,M=1.0) :)
    #Overlap_Bool=Source_Overlap_Calc(GC_RA_Arcmin,Source_RA,GC_DEC_Arcmin,Source_DEC_Arcmin,D25_Maj,0,D25_Min,0,Beta,0,M=1.0)[-1]
    ##Overlap_Bool=Source_Overlap_Calc(GC_RA_Arcmin,Source_RA_Arcmin,GC_DEC_Arcmin,Source_DEC_Arcmin,D25_Maj,0,D25_Min,0,Beta,0,M=1.0)[-1]
    #Overlap_Bool=Source_Overlap_Calc(GC_RA_Arcmin,Source_RA_Arcmin,GC_DEC_Arcmin,Source_DEC_Arcmin,D25_Maj,0,D25_Min,0,Theta,0,M=1.0)[-1]
    #Overlap_Bool=Source_Overlap_Calc(-1.0*GC_RA_Arcmin,-1.0*Source_RA_Arcmin,GC_DEC_Arcmin,Source_DEC_Arcmin,D25_Maj,0,D25_Min,0,Theta,0,M=1.0)[-1]
    Overlap_Bool=Source_Overlap_Calc(Source_RA_Arcmin,GC_RA_Arcmin,GC_DEC_Arcmin,Source_DEC_Arcmin,D25_Maj,0,D25_Min,0,Theta,0,M=1.0)[-1]
    #Overlap_Bool=Source_Overlap_Calc(Source_RA,GC_RA,GC_DEC,Source_DEC,D25_Maj,0,D25_Min,0,Theta,0,M=1.0)[-1]

    if(Overlap_Bool==True):
        return False
    if(Overlap_Bool==False):
        return True

def Circular_D25_Bool_Calc(ObsID):
    D25_Info=D25_Query(ObsID)
    D25_Maj=D25_Info[0]
    if(D25_Maj==np.nan):
        return
    #print(isinstance(D25_Query(7087)[2],np.ma.core.MaskedConstant))
    D25_Min=D25_Info[1]
    D25_Angle=D25_Info[2]
    if(isinstance(D25_Min,np.ma.core.MaskedConstant) or isinstance(D25_Angle,np.ma.core.MaskedConstant)): #In the event that the ellipse is assumed to be a cirlce by RC3
        print(str(ObsID)+" has a circular D25 region")
        return True
    else:
        return False

def Aimpoint_Coords_Calc(ObsID):
    Cur_Evt2_Filepath=File_Query(ObsID)
    hdulist = fits.open(Cur_Evt2_Filepath)
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
    Cur_Evt2_Filepath=File_Query(ObsID)
    hdulist = fits.open(Cur_Evt2_Filepath)
    #Exposure_Time=hdulist[1].header['EXPOSURE']
    Pointing_RA=hdulist[1].header['RA_PNT']
    Pointing_Dec=hdulist[1].header['DEC_PNT']
    Pointing_Diff_RA=Pointing_RA-GC_RA
    Pointing_Diff_Dec=Pointing_Dec-GC_DEC
    x1_Rad=Deg_to_Rad(GC_RA)
    x2_Rad=Deg_to_Rad(Pointing_RA)
    y1_Rad=Deg_to_Rad(GC_DEC)
    y2_Rad=Deg_to_Rad(Pointing_Dec)
    Have_Dist=Haversine_Distance(x1_Rad,x2_Rad,y1_Rad,y2_Rad)
    Have_Dist=np.abs(Have_Dist)
    Have_Dist_Deg=Rad_to_Deg(Have_Dist)
    Have_Dist_Arcmin=Have_Dist_Deg*60.0
    print("Have_Dist_Arcmin: ", Have_Dist_Arcmin)
    Dist=np.sqrt((Pointing_Diff_RA**2.0)+(Pointing_Diff_Dec**2.0))
    Dist_Arcmin=Dist*60.0
    #Dist_L=[Cur_Evt2_ObsID,Dist_Arcmin]
    #Dist_H_L.append(Dist_L)
    #return Dist_H_L
    #return Dist_Arcmin
    return Have_Dist_Arcmin
def Overlapping_ObsID_Calc(Data,Dist_Threshold=2.0):
    #pass
    #Aimpoint_Coords=np.vectorize(Aimpoint_Coords_Calc)(ObsID)
    RA_Aimpoint_A=Data["RA_Aimpoint"]
    RA_Aimpoint_L=list(RA_Aimpoint_A)
    DEC_Aimpoint_A=Data["DEC_Aimpoint"]
    DEC_Aimpoint_L=list(DEC_Aimpoint_A)
    #ObsID_A=Data['OBSID']
    ObsID_A=Data['ObsID']
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
            x1_Rad=Deg_to_Rad(RA_Aimpoint)
            x2_Rad=Deg_to_Rad(RA_Aimpoint_Test)
            y1_Rad=Deg_to_Rad(DEC_Aimpoint)
            y2_Rad=Deg_to_Rad(DEC_Aimpoint_Test)
            Have_Dist=Haversine_Distance(x1_Rad,x2_Rad,y1_Rad,y2_Rad)
            Have_Dist=np.abs(Have_Dist)
            Have_Dist_Deg=Rad_to_Deg(Have_Dist)
            Dist=Have_Dist_Deg
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
def Close_Observations_Finder(Data,ObsID,Source_Num):
    #ObsID_A=Data['OBSID']
    ObsID_A=Data['ObsID']
    ObsID_L=list(ObsID_A)
    #Source_Num_A=Data['SOURCE']
    Source_Num_A=Data['Source_Num']
    Source_Num_L=list(Source_Num_A)
    Close_ObsIDs_A=Data['Close_ObsIDs']
    Close_ObsIDs_HL=list(Close_ObsIDs_A)
    for i in range(0,len(ObsID_L)):
        ObsID_Test=ObsID_L[i]
        Source_Num_Test=Source_Num_L[i]
        if((ObsID_Test==ObsID) and (Source_Num_Test==Source_Num)):
            Match_Inx=i
            break
    Close_ObsIDs_L_Str=Close_ObsIDs_HL[i]
    Close_ObsIDs_L_Strings=re.split("[\[\];]",Close_ObsIDs_L_Str) #Note: I am not sure if this is the correct regex
    print("Close_ObsIDs_L_Strings: ", Close_ObsIDs_L_Strings)
    Close_ObsIDs_L_Strings.remove("")
    Close_ObsIDs_L=[]
    for ObsID_Str in Close_ObsIDs_L_Strings:
        if(ObsID_Str==""):
            continue
        ObsID_Int=int(ObsID_Str)
        Close_ObsIDs_L.append(ObsID_Int)
    return Close_ObsIDs_L
    #df.loc[df['column_name'].isin(some_values)]
    """
    Data_Test=Data.loc[Data['OBSID'].isin(Close_ObsIDs_L)]
    return Data_Test
    """
def Duplicate_Source_Calc(Data,Dist_Threshold=(2.0/3600.0)): #Threshold needs to be determed. Temperary value used.
    #pass
    #Aimpoint_Coords=np.vectorize(Aimpoint_Coords_Calc)(ObsID)
    #ObsID_A=Data['OBSID']
    ObsID_A=Data['ObsID']
    ObsID_L=list(ObsID_A)
    #Source_Num_A=Data['SOURCE']
    Source_Num_A=Data['Source_Num']
    Source_Num_L=list(Source_Num_A)
    Source_RA_A=Data['RA_WavD']
    Source_RA_L=list(Source_RA_A)
    Source_Dec_A=Data['DEC_WavD']
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
        print("Close_ObsIDs_L_Strings: ", Close_ObsIDs_L_Strings)
        Close_ObsIDs_L=[]
        for ObsID_Str in Close_ObsIDs_L_Strings:
            if(ObsID_Str==""):
                continue
            ObsID_Int=int(ObsID_Str)
            Close_ObsIDs_L.append(ObsID_Int)
        #df.loc[df['column_name'].isin(some_values)]
        Data_Test=Data.loc[Data['ObsID'].isin(Close_ObsIDs_L)]
        ObsID_A_Test=Data_Test['ObsID']
        ObsID_L_Test=list(ObsID_A_Test)
        Source_Num_A_Test=Data_Test['Source_Num']
        Source_Num_L_Test=list(Source_Num_A_Test)
        Source_RA_A_Test=Data_Test['RA_WavD']
        Source_RA_L_Test=list(Source_RA_A_Test)
        Source_Dec_A_Test=Data['DEC_WavD']
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
            x1_Rad=Deg_to_Rad(Source_RA)
            x2_Rad=Deg_to_Rad(Source_RA_Test)
            y1_Rad=Deg_to_Rad(Source_Dec)
            y2_Rad=Deg_to_Rad(Source_Dec_Test)
            Have_Dist=Haversine_Distance(x1_Rad,x2_Rad,y1_Rad,y2_Rad)
            Have_Dist=np.abs(Have_Dist)
            Have_Dist_Deg=Rad_to_Deg(Have_Dist)
            Dist=Have_Dist_Deg
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
def Color_Color_Calc(S,M,H):
    if(H+M==0.0):
        HC_Ratio=np.nan
    else:
        HC_Ratio=(H-M)/((H+M))
    if(M+S==0.0):
        SC_Ratio=np.nan
    else:
        SC_Ratio=(M-S)/((M+S))
    return HC_Ratio, SC_Ratio
def Flux_Colors_Calc(ObsID, Source_Num, Data_Path="/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All.csv", Soft_Key="NET_FLUX_APER_0.3-1.0", Medium_Key="NET_FLUX_APER_1.0-2.1", Hard_Key="NET_FLUX_APER_2.1-7.5"):
    Glob_L=glob.glob(Data_Path)
    if(len(Glob_L)!=1):
        raise Exception("Data_Path has "+str(len(Glob_L))+" Matching Files ! ! !")
    Data_Path=Glob_L[0]
    Data=pd.read_csv(Data_Path)
    Data_Reduced=Data[Data["ObsID"] == ObsID]
    Data_Reduced=Data_Reduced.reset_index(drop=True)
    Soft_A=Data_Reduced[Soft_Key]
    Medium_A=Data_Reduced[Medium_Key]
    Hard_A=Data_Reduced[Hard_Key]
    Soft=Soft_A.iloc[Source_Num-1]
    Medium=Medium_A.iloc[Source_Num-1]
    Hard=Hard_A.iloc[Source_Num-1]
    print("Soft: ", Soft)
    print("Medium: ", Medium)
    print("Hard: ", Hard)
    Colors=Color_Color_Calc(Soft,Medium,Hard)
    return Colors
def CCD_Region_Query(ObsID):
    #Region_File_Query_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+str(Gname)+"/Area_Lists/"+str(ObsID)+"/*"+str(ObsID)+"*_CCD_Regions_simple_region_modifed_Code.txt"
    Region_File_Query_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/*/Area_Lists/"+str(ObsID)+"/*"+str(ObsID)+"*_CCD_Regions_simple_region_modifed_Code.txt"
    Region_File_Path_L=glob.glob(Region_File_Query_Path)
    if(len(Region_File_Path_L)!=1):
        #print str(ObsID)+" has "+str(len(evtfpath_L)-1)+"Additional Evt2 Files ! ! !"
        print(str(ObsID)+" has "+str(len(Region_File_Path_L))+" Simple Region Files ! ! !")
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
        #print "Point: ", Point
        X=Point[0]
        Y=Point[1]
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
def Region_Lengths_Converter(Lengths_L,To_RA_Dec_Bool=False):
    #2.03252032520325 the converstion factor is 2.03252032520325pix/arcsec (7317.0731707317 pix/deg) or 0.4920+-0.0001 arcsec/pix
    Lengths_A=np.array(Lengths_L)
    if(To_RA_Dec_Bool):
        Lengths_Coverted_A=(0.4920*Lengths_A)/3600.0
    else:
        Lengths_Coverted_A=7317.0731707317*Lengths_A #7317.0731707317 pix/deg
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
        #return CXCRegion(Region_L[0])
    Unioned_Region=Region_L[0]
    #Unioned_Region=CXCRegion(Unioned_Region)
    print("Initial Unioned_Region: ", Unioned_Region)
    print("type(Unioned_Region): ", type(Unioned_Region))
    for i in range(1,len(Region_L)):
        Region=Region_L[i]
        #Region=CXCRegion(Region)
        #print "Loop Region: ", Region
        #print "type(Region): ", type(Region)
        Unioned_Region=Unioned_Region+Region
    return Unioned_Region
def Region_Intersection(Region_L):
    if(len(Region_L)==0):
        return
    if(len(Region_L)==1):
        return Region_L[0]
        #return CXCRegion(Region_L[0])
    Intersected_Region=Region_L[0]
    #print "str(Intersected_Region.shapes[0]): ", str(Intersected_Region.shapes)
    #Intersected_Region=CXCRegion(str(Intersected_Region.shapes))
    for i in range(1,len(Region_L)):
        Region=Region_L[i]
        #Region=CXCRegion(Region)
        Intersected_Region=Intersected_Region*Region
    return Intersected_Region
def Limiting_Flux_Annulus_Intersection(ObsID,Source_Num):
    #pass
    CCD_Region_Fpath=CCD_Region_Query(ObsID)
    X_Aimpoint,Y_Aimpoint=Aimpoint_Coords_Calc(ObsID)
    Annulus_Number=Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num)
    Annulus_Shape_Str=Annulus_Shape_String_Generator(X_Aimpoint,Y_Aimpoint,Annulus_Number,W=(1.0/60.0))
    Annulus_Shape_Data=Region_Parse(Annulus_Shape_Str)[0]
    #print "Annulus_Shape_Data: ", Annulus_Shape_Data
    Annulus_Shape_Point_L=Annulus_Shape_Data[1]
    #print "Annulus_Shape_Point_L: ", Annulus_Shape_Point_L
    Annulus_Shape_Radii_L=Annulus_Shape_Data[2]
    #print "Annulus_Shape_Radii_L: ", Annulus_Shape_Radii_L
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
def Region_ObsID_Reference_Frame_Converter(ObsID,Region):  #Note: Still being worked on
    #pass
    Region_Str=''
    print("Region: ", Region)
    Region_Shapes_L=Region.shapes
    #Reg_Combined=""
    #for Reg in Region:
    #for Reg in Region_Shapes_L:
    for i in range(0,len(Region_Shapes_L)):
        Reg=Region_Shapes_L[i]
        print("Reg: ", Reg)
        print("type(Reg):", type(Reg))
        #Reg_Str=Reg.shapes
        #print "Reg_Str: ", Reg_Str
        Reg_Str=str(Reg)
        print("Reg_Str: ", Reg_Str)
        if(i==0):
            Reg_Combined=Region_Shapes_L[0]
        """
        if(Reg_Str=='annulus')
            Reg_Parse_Data=Region_Parse(Reg_Str)[0]
            print "Reg_Parse_Data: ", Reg_Parse_Data
            Reg_Parse_Shape=Reg_Parse_Data[0]
            print "Reg_Parse_Shape: ", Reg_Parse_Shape
            Reg_Parse_Point_L=Reg_Parse_Data[1]
            print "Reg_Parse_Point_L: ", Reg_Parse_Point_L
            Reg_Parse_Length_L=Reg_Parse_Data[2]
            print "Reg_Parse_Length_L: ", Reg_Parse_Length_L
        """
        X_L=Reg.xpoints
        Y_L=Reg.ypoints
        #print "Y_L: ", Y_L
        Length_L=Reg.radii
        Angles_L=Reg.angles
        Logic_Operator=Reg.logic
        print("Logic_Operator: ", Logic_Operator)
        print("type(Logic_Operator): ", type(Logic_Operator))
        Logic_Operator_Str=str(Logic_Operator.str)
        print("Logic_Operator_Str: ", Logic_Operator_Str)
        Point_L=[]
        for i in range(0,len(X_L)):
            Cur_Point=[X_L[i],Y_L[i]]
            Point_L.append(Cur_Point)
        Point_L_Converted=Region_Coords_Converter(ObsID,Point_L)
        Length_L_Converted=Region_Lengths_Converter(Length_L)
        Point_Converted=Point_L_Converted[0]
        X_Converted=Point_Converted[0]
        Y_Converted=Point_Converted[1]
        #if(Reg=='annulus'):
        if('annulus' in Reg_Str):
            Inner_Radius_Converted=Length_L_Converted[0]
            Outer_Radius_Converted=Length_L_Converted[1]
            Str_Converted='annulus(' + str(X_Converted) +','+ str(Y_Converted)+','+ str(Inner_Radius_Converted)+','+ str(Outer_Radius_Converted)+')'
        #if(('box' in Reg) or (Reg=='rotbox')):
        if('box' in Reg_Str):
            Width=Length_L_Converted[0]
            Hight=Length_L_Converted[1]
            Angle=Angles_L[0] #Note: This might be incorrect. The two differnt observations can have two differnt roll angles and thus would require the angle to be converted. I have to make sure this is the case though.
            Str_Converted='box(' + str(X_Converted) +','+ str(Y_Converted)+','+ str(Width)+','+ str(Hight)+','+str(Angle)+')'
        #Region_Str=Region_Str+Str_Converted+'\n' #Note: This is inccorect. The combined regions must be constructed with logic operators (as intersection ("*") or union ("+"))
        #Region_Str=Region_Str Logic_Operator Str_Converted
        Region_Str=Region_Str+Logic_Operator_Str+Str_Converted #Note: This does not work. It might be better to convert into a region objects frist and then combine the regions objects rather than trying to convert the combined string region.
        """
        #print "Region_Str: ", Region_Str
        print "Str_Converted:", Str_Converted
        Cur_Region=CXCRegion(Str_Converted)
        if(i!=0):
            if(Logic_Operator_Str=="+"):
                Reg_Combined=Reg_Combined+Cur_Region
            if(Logic_Operator_Str=="*"):
                Reg_Combined=Reg_Combined*Cur_Region
        """
    print("Region_Str Final: ", Region_Str)
    Region_Converted=CXCRegion(Region_Str) #Note: Should use CXCRegion objects convenience functions instead
    #Region_Converted=Reg_Combined
    return Region_Converted
def Unioned_Limiting_Flux_Area_Calc(Data,F_L_Bounds): #Note: Still being worked on
    #pass
    ObsID_A=Data['ObsID']
    Src_A=Data['Source_Num']
    Data["Annulus_Intersection_Regions"]=np.vectorize(Limiting_Flux_Annulus_Intersection)(ObsID_A,Src_A)
    print('Data["Annulus_Intersection_Regions"]: ', Data["Annulus_Intersection_Regions"])
    ##Data_Slice=Limiting_Flux_Data_Slice_Calc(Data,F_L_Bounds)
    Data_Slice=Data #Note: This is a test
    Data_Slice=Data_Slice.drop_duplicates(subset=["ObsID","Offaxis_Angle_Annulus_Number"])
    Data_Slice=Data_Slice[Data_Slice['Close_ObsIDs_Bool']]
    Data_Slice=Data_Slice.reset_index(drop=True)
    print("Data_Slice:\n", Data_Slice)
    ObsID_A_Slice=Data_Slice['ObsID']
    ObsID_L_Slice=list(ObsID_A_Slice)
    Src_A_Slice=Data_Slice['Source_Num']
    Src_L_Slice=list(Src_A_Slice)
    for i in range(0,len(ObsID_L_Slice)):
        ObsID=ObsID_L_Slice[i]
        Source_Num=Src_L_Slice[i]
        if(ObsID not in list(Data_Slice['ObsID'])):
            pass
        Close_ObsIDs_L=Close_Observations_Finder(Data_Slice,ObsID,Source_Num)
        print("Close_ObsIDs_L: ", Close_ObsIDs_L)
        #Example: Data_Test=Data.loc[Data['OBSID'].isin(Close_ObsIDs_L)]
        Data_Slice_Close=Data_Slice.loc[Data_Slice['ObsID'].isin(Close_ObsIDs_L)] #Note: There is a bug here
        print("Data_Slice_Close:\n", Data_Slice_Close)
        ObsID_L_Slice_Close=list(Data_Slice_Close['ObsID'])
        Data_Slice=Data_Slice.loc[~Data_Slice['ObsID'].isin(Close_ObsIDs_L)] #Note: This is to prevent the repeated calculation of the area of the same Overlapping ObsID Group
        Annulus_Intersection_Regions_A=Data_Slice_Close["Annulus_Intersection_Regions"] #Note: All Annulus_Intersection_Regions must be put in the same ObsID reference frame
        Annulus_Intersection_Regions_L=list(Annulus_Intersection_Regions_A)
        Annulus_Intersection_Regions_Converted_L=[]
        for j in range(0,len(ObsID_L_Slice_Close)):
            ObsID_Test=ObsID_L_Slice_Close[j]
            Annulus_Intersection_Region=Annulus_Intersection_Regions_L[j]
            if(ObsID_Test==ObsID):
                Annulus_Intersection_Regions_Converted_L.append(Annulus_Intersection_Region)
                continue
            #Annulus_Intersection_Region_Converted=Region_ObsID_Reference_Frame_Converter(ObsID_Test,Annulus_Intersection_Region)
            Annulus_Intersection_Region_Converted=Region_ObsID_Reference_Frame_Converter(ObsID,Annulus_Intersection_Region)
            Annulus_Intersection_Regions_Converted_L.append(Annulus_Intersection_Region_Converted)
        print("Annulus_Intersection_Regions_Converted_L: ", Annulus_Intersection_Regions_Converted_L)
        Unioned_Region=Region_Union(Annulus_Intersection_Regions_Converted_L)
        print("Unioned_Region: ",Unioned_Region)
        Unioned_Area=Unioned_Region.area()
        return Unioned_Area
def Source_Counts_To_Flux_Converter(fpath,Outfpath,CR_K=0.01):
    data=pd.read_csv(fpath)
    Key_L=list(data.keys())
    print("Key_L: ", Key_L)
    """
    data['RAW_COUNTS[.3-7.5]']=data['RAW_COUNTS[.3-1]']+data['RAW_COUNTS[1-2.1]']+data['RAW_COUNTS[2.1-7.5]']
    data['BKG_COUNTS[.3-7.5]']=data['BKG_COUNTS[.3-1]']+data['BKG_COUNTS[1-2.1]']+data['BKG_COUNTS[2.1-7.5]']
    data['NET_COUNTS[.3-7.5]']=data['RAW_COUNTS[.3-7.5]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[.3-7.5]'])
    data['NET_COUNTS[.3-1]']=data['RAW_COUNTS[.3-1]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[.3-1]'])
    data['NET_COUNTS[1-2.1]']=data['RAW_COUNTS[1-2.1]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[1-2.1]'])
    data['NET_COUNTS[2.1-7.5]']=data['RAW_COUNTS[2.1-7.5]']-((data['AREA']/data['BKG_AREA'])*data['BKG_COUNTS[2.1-7.5]'])
    """
    #Luminosity_Calc
    #ObsID_A=data['OBSID']
    ObsID_A=data['ObsID']
    ##'''
    Exposure_Time_A=ObsID_A.apply(Exposure_Time_Calc) #This works but takes to long as it unnecessarily repeats getting the exposure time for every object rather than every ObsID
    #print "Exposure_Time_A:\n", Exposure_Time_A
    data['Exposure_Time']=Exposure_Time_A
    #Start_Date_Calc
    Start_Date_A=ObsID_A.apply(Start_Date_Calc)
    data['Start_Date']=Start_Date_A
    #data['Count_Rate[.3-7.5]']=data['RAW_COUNTS[.3-7.5]']/data['Exposure_Time'] #This should use net counts instead of raw counts!
    #data['Count_Rate[.3-7.5]']=data['NET_COUNTS[.3-7.5]']/data['Exposure_Time']
    #Src_A=data['SOURCE']
    Src_A=data['Source_Num']
    #df['new_column'] = np.vectorize(fx)(df['A'], df['B'])
    #data['Known_Flux[.3-7.5]']=np.vectorize(Source_Known_Flux_Finder,excluded=['Effective_Area_Correction_Bool'])(ObsID_A,Src_A,Effective_Area_Correction_Bool=False) #Note: This is done because the counts from the counts_info.csv file are already effective area corrected. I am not sure if this is the correct method since this also effects the pile-up estimation. Therefore this is an improvement of the double effective area correction bug but not a full solution. The fluxes will still be incorrect!
    #print "data:\n", data
    """
    print "data Grouped:\n", data.groupby("OBSID").groups
    Data_Grouped=data.groupby("OBSID")
    for name, group in Data_Grouped:
        print(name)
        print(group)
    """
    #F_U=F_K*((float(C)/Exposure_Time)/CR_K) # F_U:-float, Flux Unknown, This converts the current count value to the current flux value by assuming a linear relationship bewteen counts and flux, F_U is the calculated flux
    #data['Flux[.3-7.5]']=data['Known_Flux[.3-7.5]']*(data['Count_Rate[.3-7.5]']/CR_K)
    ##data['Limiting_Flux[.3-7.5]']=np.vectorize(Limiting_Flux_Calc)(ObsID_A,Src_A)
    #data['Flux_Cut_Bool[.3-7.5]']=np.where(data['Flux[.3-7.5]']>data['Flux_Cut_Bool[.3-7.5]'])
    #Flux_A=data['Flux[.3-7.5]']
    #print "type(Flux_A[3]): ", type(Flux_A[3])
    ##Limiting_Flux_A=data['Limiting_Flux[.3-7.5]']
    #print "type(Limiting_Flux_A[3]): ", type(Limiting_Flux_A[3])
    ##data['Flux_Cut_Bool[.3-7.5]']=np.vectorize(Flux_Cut_Bool_Calc)(Flux_A,Limiting_Flux_A)
    data["Inside_FOV_Bool"]=np.vectorize(Inside_FOV_Bool_Calc)(ObsID_A,Src_A)
    Coords=np.vectorize(RA_DEC_Calc)(ObsID_A,Src_A)
    data["RA_WavD"]=Coords[0]
    data["DEC_WavD"]=Coords[1]
    #Offaxis_Angle_Calc(ObsID,Source_Num)
    data["Offaxis_Angle"]=np.vectorize(Offaxis_Angle_Calc)(ObsID_A,Src_A)
    #Offaxis_Angle_Annulus_Number_Calc(ObsID,Source_Num)
    data["Offaxis_Angle_Annulus_Number"]=np.vectorize(Offaxis_Angle_Annulus_Number_Calc)(ObsID_A,Src_A)
    ##data["Offaxis_Angle_Annulus_CCD_Incompleteness"]=np.vectorize(Offaxis_Angle_Annulus_CCD_Incompleteness_Calc)(ObsID_A,Src_A)
    data["Offaxis_Angle_Annulus_Area"]=np.vectorize(Offaxis_Angle_Annulus_Area_Calc)(ObsID_A,Src_A)
    ##data["Observed_Offaxis_Angle_Annulus_Area"]=np.vectorize(Observed_Offaxis_Angle_Annulus_Area_Calc)(ObsID_A,Src_A)
    data["Gname_Modifed"]=np.vectorize(Gname_Query)(ObsID_A)
    #Galaxy_Morph_Query(ObsID)
    data["Galaxy_Morph"]=np.vectorize(Galaxy_Morph_Query)(ObsID_A)
    data["Galaxy_Morph_Simple"]=np.vectorize(Morph_Reducer)(ObsID_A)
    #GC_Query(ObsID)
    GC_Coords=np.vectorize(GC_Query)(ObsID_A)
    data["RA_GC"]=GC_Coords[0]
    data["DEC_GC"]=GC_Coords[1]
    #Galatic_Distance_Query
    data["Galatic_Distance"]=np.vectorize(Galatic_Distance_Query)(ObsID_A)
    #D25_Query(ObsID)
    ##data["D25"]=np.vectorize(D25_Query)(ObsID_A)
    D25_Tuple_A=np.vectorize(D25_Query)(ObsID_A)
    data["D25"]=D25_Tuple_A[0]
    data["D25_Maj"]=D25_Tuple_A[0]
    data["D25_Min"]=D25_Tuple_A[1]
    data["D25_Angle"]=D25_Tuple_A[2]
    #Distance_From_GC_Calc(ObsID,Source_Num)
    data["Source_Distance_From_GC"]=np.vectorize(Distance_From_GC_Calc)(ObsID_A,Src_A)
    #Psi_Calc
    data["Psi"]=np.vectorize(Psi_Calc)(ObsID_A,Src_A)
    #Outside_D25_Bool_Calc(ObsID,Source_Num)
    data["Outside_D25_Bool"]=np.vectorize(Outside_D25_Bool_Calc)(ObsID_A,Src_A)
    #Outside_Elliptical_D25_Bool_Calc
    data["Outside_Elliptical_D25_Bool"]=np.vectorize(Outside_Elliptical_D25_Bool_Calc)(ObsID_A,Src_A)
    #Circular_D25_Bool_Calc(ObsID)
    data["Circular_D25_Bool"]=np.vectorize(Circular_D25_Bool_Calc)(ObsID_A)
    #Distance_Galatic_Center_to_Aimpoint_Calc(ObsID)
    ##'''
    #Aimpoint_Coords_Calc(ObsID)
    Aimpoint_Coords=np.vectorize(Aimpoint_Coords_Calc)(ObsID_A)
    data["RA_Aimpoint"]=Aimpoint_Coords[0]
    data["DEC_Aimpoint"]=Aimpoint_Coords[1]
    data["Roll_Angle"]=np.vectorize(Roll_Angle_Calc)(ObsID_A)
    data["Distance_GC_to_Aimpoint"]=np.vectorize(Distance_Galatic_Center_to_Aimpoint_Calc)(ObsID_A)
    data=Luminosity_Calc_Bulk(data)
    data=Overlapping_ObsID_Calc(data)
    data=Duplicate_Source_Calc(data)
    #with pd.option_context('display.max_rows', None):  # more options can be specified also
        #print "data:\n", data
    Header_L = data.columns.tolist()
    print("Header_L: ", Header_L)
    #Header_Modified_L=Header_L[:21]+Header_L[281:]+Header_L[21:281]
    #data=data[Header_Modified_L]
    print("data:\n", data)
    #Outfname=fpath.split("/")[0].split(".")[0]+"_Flux_Calc.csv"
    #print "Outfname: ", Outfname
    #Outpath=fpath.split("/")[:-1]

    data.to_csv(Outfpath)
    #print "data.dtypes:\n", data.dtypes
    #print "data:\n", data.where(data['Flux_Cut_Bool[.3-7.5]']==True,True,False)
    #print "data:\n", data.loc(data['Flux_Cut_Bool[.3-7.5]']==True)
    #print "data:\n", data.loc(True)
    #print data['Flux_Cut_Bool[.3-7.5]'].loc(True)
def Test_Column_Mods(fpath="/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All_Modified.csv"):
    data=pd.read_csv(fpath)
    #print("data: ", data)
    Header_L = data.columns.tolist()
    #print("Header_L: ", Header_L)
    #Header_Modified_L=Header_L[:21]+Header_L[281:]+Header_L[21:281]
    print("Header_L[:21]: ", Header_L[:21])
    print("Header_L[282:]: ", Header_L[282:])
    Header_Modified_L=Header_L[:21]+Header_L[282:]+Header_L[21:282]
    #print("Header_Modified_L: ", Header_Modified_L)
    data=data[Header_Modified_L]
    Header_L_Modified = data.columns.tolist()
    ##print("Header_L_Modified: ", Header_L_Modified)
def Data_Testing(Func,fpath,arg,Outside_D25_Bool=False):
    data=pd.read_csv(fpath)
    if(Outside_D25_Bool==True):
        #data=data.drop(data[data.Outside_D25_Bool==False].index, inplace=True)
        #new_dataframe = a_dataframe[a_dataframe.B <= 3]
        data=data[data.Outside_D25_Bool]
        data=data.reset_index(drop=True)
    Output=Func(data,arg) #Note: This is a temperay fix for *args
    return Output
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
    print("Slope_A: ",Slope_A)
    ##Flux_Error_A=np.linspace(Flux_Error_Bounds[0],Flux_Error_Bounds[1],Steps)
    Flux_Error_A=np.geomspace(Flux_Error_Bounds[0],Flux_Error_Bounds[1],Steps)
    print("Flux_Error_A: ", Flux_Error_A)
    F_L_Bounds_Full=[5e-17,2e-15]
    ##F_L_Bounds_A=np.linspace(F_L_Bounds_Full[0],F_L_Bounds_Full[1],10)
    ##F_L_Bounds_A=np.linspace(F_L_Bounds_Full[0],F_L_Bounds_Full[1],20)
    F_L_Bounds_A=np.geomspace(F_L_Bounds_Full[0],F_L_Bounds_Full[1],20)
    F_L_Bounds_HL=[]
    for i in range(len(F_L_Bounds_A)-1):
        Cur_Bounds=[F_L_Bounds_A[i],F_L_Bounds_A[i+1]]
        F_L_Bounds_HL.append(Cur_Bounds)
    print("F_L_Bounds_HL: ", F_L_Bounds_HL)
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
    print("Cur_Best_Flux_Cut_Frac_L: ", Cur_Best_Flux_Cut_Frac_L)
    return [Cur_Best_Fit_Parameters,Min_Modified_Standard_Dev]

def Limiting_Flux_Area_Calc(Data,F_L_Bounds):
    Data_Slice=Limiting_Flux_Data_Slice_Calc(Data,F_L_Bounds)
    Data_Slice=Data_Slice.drop_duplicates(subset=["ObsID","Offaxis_Angle_Annulus_Number"])
    Data_Slice=Data_Slice.reset_index(drop=True)
    print("F_L_Bounds: ", F_L_Bounds)
    #with pd.option_context('display.max_columns', None):
    #print "Limiting_Flux_Area Data_Slice:\n", Data_Slice
    print("Limiting_Flux_Area Data_Slice:\n", Data_Slice[["ObsID","Gname_Modifed","RA","Offaxis_Angle_Annulus_Number"]])
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
    print("data.columns:\n", data.columns)
    #print "data:\n", data
    with pd.option_context('display.max_columns', None):  # more options can be specified also
        print("data:\n", data)
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data)
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Slope_Bounds=[0,100.0],Flux_Error_Bounds=[0,5])
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Flux_Error_Bounds=[1e-16,20e-16])
    ##Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Flux_Error_Bounds=[1e-16,40e-16])
    if(Flux_Model_B):
        Flux_Cut_Fit=Limiting_Flux_Model_Fitting(data,Flux_Error_Bounds=[1e-16,1e-13])
        print("Flux_Cut_Fit: ", Flux_Cut_Fit)
    Test_Bool=False
    Flux_A=data['Flux[.3-7.5]']
    #Flux_A=Flux_A/10.0 #Note: This is a test. THIS CREATES INCORRECT DATA!!!! DO NOT USE FOR SCIENCE!!!!
    #Flux_A=Flux_A/5.0 #Note: This is a test. THIS CREATES INCORRECT DATA!!!! DO NOT USE FOR SCIENCE!!!!
    Flux_Avg=Flux_A.mean()
    print("Flux_Avg: ", Flux_Avg)
    Flux_L=list(Flux_A)
    Limiting_Flux_A=data['Limiting_Flux[.3-7.5]']
    Limiting_Flux_Avg=Limiting_Flux_A.mean()
    print("Limiting_Flux_Avg: ", Limiting_Flux_Avg)
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
    print(Test_Bool)
    print(n)
    if(Flux_Model_B):
        Model_Limiting_Flux_A=Limiting_Flux_Model(Limiting_Flux_A,Flux_Cut_Fit[0][0],Flux_Cut_Fit[0][1])
        print("Model_Limiting_Flux_A: ", Model_Limiting_Flux_A)
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
    print("LogN_LogS_Hist:\n", LogN_LogS_Hist)
    LogN_LogS_Bins=LogN_LogS_Hist_A[1]
    print("LogN_LogS_Bins:\n", LogN_LogS_Bins)
    print("type(LogN_LogS_Bins): ", type(LogN_LogS_Bins))
    Limiting_Flux_Area_L=[]
    for i in range(0,len(list(LogN_LogS_Bins))-1):
        j=i+1
        F_L_Bin_Bounds=[LogN_LogS_Bins[i],LogN_LogS_Bins[j]]
        Limiting_Flux_Area=Limiting_Flux_Area_Calc(data,F_L_Bin_Bounds)
        #print "Limiting_Flux_Area: ", Limiting_Flux_Area
        Limiting_Flux_Area_L.append(Limiting_Flux_Area)
    Limiting_Flux_Area_A=np.array(Limiting_Flux_Area_L)
    print("Limiting_Flux_Area_A:\n", Limiting_Flux_Area_A)
    Limiting_Flux_Area_A_Summed=np.cumsum(Limiting_Flux_Area_A[::-1])[::-1]
    print("Limiting_Flux_Area_A_Summed:\n", Limiting_Flux_Area_A_Summed)
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
    print("LogN_LogS_Source_Density_Hist_A:\n", LogN_LogS_Source_Density_Hist_A)
    LogN_LogS_Bins_A=np.array(LogN_LogS_Bins)
    LogN_LogS_Bins_Left_Edges_A=LogN_LogS_Bins_A[:-1]
    #plt.bar(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A)
    ##plt.plot(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A,label="Data")
    plt.semilogx(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A,label="Data",marker=".")
    ##plt.errorbar(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Source_Density_Hist_A,yerr=LogN_LogS_Source_Density_Error_A,linestyle="") #Current Functioning Error bars
    ##plt.semilogx(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Hist/5.0,label="Test",marker=".")
    LogN_LogS_Soft_A=Background_Source_Calc(LogN_LogS_Bins_Left_Edges_A,Hardness_Str="S")
    #plt.plot(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Soft_A,label="Soft_Giacconi",linestyle="--")
    plt.semilogx(LogN_LogS_Bins_Left_Edges_A,LogN_LogS_Soft_A,label="Soft_Giacconi",linestyle="--")
    LogN_LogS_Hard_A=Background_Source_Calc(LogN_LogS_Bins_Left_Edges_A,Hardness_Str="H")
    print("LogN_LogS_Hard_A:\n", LogN_LogS_Hard_A)
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
#print(Gname_Query(10125))
#print GC_Query(12095)
#print GC_Query(6869)
#5905
#print GC_Query(5905)
#9552
#print GC_Query(9552)
#print D25_Query(12095)
#print Source_Known_Flux_Finder(12095,1)
#print Source_Known_Flux_Finder(12095,1,False)
#print Limiting_Flux_Calc(12095,1)
#print Limiting_Flux_Calc(12095,5)
#print Limiting_Flux_Calc(2197,5)
#print Limiting_Flux_Calc(6869,1)
#print(Distance_From_GC_Calc(12095,5))
#print(Distance_From_GC_Calc(735,7))
#print Distance_From_GC_Calc(6869,5)
#print Outside_D25_Bool_Calc(12095,5)
#print Distance_From_GC_Calc(12095,20)
#print Outside_D25_Bool_Calc(12095,20)
#print(Distance_Galatic_Center_to_Aimpoint_Calc(12095))
#print(Distance_Galatic_Center_to_Aimpoint_Calc(6869))
#print(Distance_Galatic_Center_to_Aimpoint_Calc(735))
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
##Source_Counts_To_Flux_Converter("/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info.csv","/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv")
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv')
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv',Outside_D25_Bool=True)
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv')
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True)
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True,Output_File_Ext="png")
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_Flux_Calc.csv',Outside_D25_Bool=True)
#Data_Analysis('/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/Old_Runs/Old_092720/counts_info_Flux_Calc.csv',Outside_D25_Bool=True)
#print Data_Testing(Unioned_Limiting_Flux_Area_Calc,'/Volumes/xray/anthony/Simon_Sandboxed_Code/Source_Counts_To_Flux_Converter/counts_info_testing_small_Flux_Calc.csv',[1.4E-16,1.7E-16])

#print(File_Query(10125))
#print(Gname_Query(10125))
#print(GC_Query(10125))
#print(Galaxy_Morph_Query(10125))
#print(Galaxy_Morph_Query(808))
#print(D25_Query(10125))
#10725
#print(D25_Query(10725))
#NGC 3344
#print(D25_Query(7087))
#print(type(D25_Query(7087)[2]))
#print(str(D25_Query(7087)[2]))
#print(isinstance(D25_Query(7087)[2],np.ma.core.MaskedConstant))
#print(type(D25_Query(7087)[0]))
#7850
#print(D25_Query(7850))
#print(D25_Query(10293))
#print(Exposure_Time_Calc(10125))
#print(Roll_Angle_Calc(10125))
#print(Inside_FOV_Bool_Calc(10125,1))
#print(Offaxis_Angle_Calc(10125,1))
#print(Offaxis_Angle_Annulus_CCD_Incompleteness_Calc(10125,15)) #Not working for all ObsIDs
#print(RA_DEC_Calc(10125,1))
#print(Distance_From_GC_Calc(10125,1))
#print(Outside_D25_Bool_Calc(10125,1))
#print(Aimpoint_Coords_Calc(10125))
#print(Distance_Galatic_Center_to_Aimpoint_Calc(10125))
#GC_Query(20495)
#print(D25_Query(354))
#print(D25_Query(7069))
#Source_Counts_To_Flux_Converter("/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All.csv","/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All_Modified.csv")
#Source_Counts_To_Flux_Converter("/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All.csv","/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All_Modified_2.csv")
##Source_Counts_To_Flux_Converter("/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All.csv","/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_All_Modified_3.csv")
#/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_Test.csv
#Source_Counts_To_Flux_Converter("/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_Test_2.csv","/opt/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_Test_2_Modified.csv")
#Test_Column_Mods()
#print(Phi_Query(10125, 1))
#print(Psi_Calc(10125, 1))
#print(Psi_Calc(16260, 59))
#print(Outside_Elliptical_D25_Bool_Calc(10125, 1))
#print(Outside_Elliptical_D25_Bool_Calc(10125, 15))
#print(Outside_Elliptical_D25_Bool_Calc(969, 86))
#print(Outside_Elliptical_D25_Bool_Calc(969, 87))
#print(Outside_Elliptical_D25_Bool_Calc(969, 5))
#print(Outside_Elliptical_D25_Bool_Calc(969, 88))
#print(Outside_Elliptical_D25_Bool_Calc(969, 80))
#print(Outside_Elliptical_D25_Bool_Calc(969, 84))
#print(Outside_Elliptical_D25_Bool_Calc(735, 6))
#print(Outside_Elliptical_D25_Bool_Calc(735, 7))
#Note: ObsID 17032 is on the has GC_Dec=-0.5 deg
#Note: ObsID 18760 is on the has GC_Dec=0.1 deg
#Note: ObsID 4694 is on the has GC_Dec=-0.036 deg
#print(Outside_Elliptical_D25_Bool_Calc(4694, 5))
#print(Outside_Elliptical_D25_Bool_Calc(4694, 47))
#print(Outside_Elliptical_D25_Bool_Calc(4694, 48))
#print(Outside_Elliptical_D25_Bool_Calc(4694, 46))
#print(Outside_Elliptical_D25_Bool_Calc(4694, 25))
#print(Outside_Elliptical_D25_Bool_Calc(4694, 44))
#print(Outside_Elliptical_D25_Bool_Calc(4694, 3))
#print(Galatic_Distance_Query(10125))
#print(Galatic_Distance_Query(735))
#print(Circular_D25_Bool_Calc(10125))
#print(Circular_D25_Bool_Calc(7087))
#print(Flux_Colors_Calc(10125, 1))
print(Flux_Colors_Calc(10125, 3))
