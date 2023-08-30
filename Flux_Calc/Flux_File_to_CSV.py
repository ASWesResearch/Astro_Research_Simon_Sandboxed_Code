import pandas as pd
from astropy.io import fits
from astropy.table import Table
import glob
def Flux_File_to_CSV(Path, Outpath):
    with fits.open(Path, memmap=True) as hdu_list: # set memmap=True for large files
        hdu = hdu_list[1]
        table = Table(hdu.data)
        names = [name for name in table.colnames if len(table[name].shape) <= 1]
        #df = table.to_pandas()
        #tbl[names].to_pandas(...)
        df=table[names].to_pandas()
        print("Data_Frame:\n", df)
        #print(table[0])
        #table=table[0]
        #table.write('Test_Table_3.ecsv', delimiter='\t', format='ascii')
        #table.write('Test_Table_3.ecsv', delimiter=' ', format='ascii')
        df.to_csv(Outpath)

#Flux_File_to_CSV("/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/10125/Source_Flux_10125_1_soft.flux")

def Flux_File_to_CSV_Big_Input(ObsID="All",Source_Num="All",Key="All"):
    if(Key=="All"):
        Key="*"
    if(Source_Num=="All"):
        Source_Num="*"
    if(ObsID=="All"):
        ObsID="*"
    ##Glob_Str="/opt/xray/anthony/expansion_backup/Source_Fluxes/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_"+str(Source_Num)+"_"+Key+".flux" #Real One
    #/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/10125
    Glob_Str="/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_"+str(Source_Num)+"_"+Key+".flux" #For Testing
    print("Glob_Str: ", Glob_Str)
    Path_L=glob.glob(Glob_Str)
    print("Path_L: ", Path_L)
    print("len(Path_L): ", len(Path_L))
    for Path in Path_L:
        Path_Seg_L=Path.split("_")
        print("Path_Seg_L: ", Path_Seg_L)
        Detected_Key=Path_Seg_L[-1].split(".flux")[0]
        print("Detected_Key: ", Detected_Key)
        #Detected_Source_Num=int(Path_Seg_L[-2])
        Path_No_Extension=Path.split(".flux")[0]
        #print("Path_No_Extension: ", Path_No_Extension)
        Outpath=Path_No_Extension+".csv"
        print("Outpath: ", Outpath)
        Flux_File_to_CSV(Path, Outpath)

def Flux_File_to_CSV_Big_Input_L(ObsID_L):
    for ObsID in ObsID_L:
        Flux_File_to_CSV_Big_Input(ObsID=ObsID)

def Find_Source_Num(Fpath):
    Fname=Fpath.split("/")[-1]
    Src_Num_Str=Fname.split("_")[3]
    Src_Num=int(Src_Num_Str)
    return Src_Num

def Find_Source_Num_Str(Fpath):
    Fname=Fpath.split("/")[-1]
    Src_Num_Str=Fname.split("_")[3]
    return Src_Num_Str

def Merge_CSVs(ObsID):
    #Glob_Str="/opt/xray/anthony/expansion_backup/Source_Fluxes/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_*.csv"
    #Glob_Str="/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_*.csv" #For Testing
    ###
    #DEFINE KEY LIST HERE Key_L
    ###
    Bands_Str_Raw="0.3-1,1-2,2-8,0.3-8,1-2.1,2.1-7.5,0.3-7.5"
    Bands_Str_L=Bands_Str_Raw.split(",")
    print("Bands_Str_L: ", Bands_Str_L)
    Bands_L=[]
    for Band_Str in Bands_Str_L:
        Band_Str_L=Band_Str.split("-")
        print("Band_Str_L: ", Band_Str_L)
        Band_L=[]
        for Energy_Str in Band_Str_L:
            Energy=float(Energy_Str)
            Band_L.append(Energy)
        Band_L_Med=(Band_L[0]+Band_L[1])/2.0
        #Band_L.append(Band_L_Med)
        Bands_L.append(Band_L)
    print("Bands_L: ", Bands_L)
    #"0.5:7:2.5,ultrasoft"
    Bands_Str=""
    Bands_Str_L_Reduced=[]
    for Band_L in Bands_L:
        Cur_Band_Str=str(Band_L[0])+"-"+str(Band_L[1])
        Bands_Str_L_Reduced.append(Cur_Band_Str)
    #print("Bands_Str_L_Reduced: ", Bands_Str_L_Reduced)
    CSC_Keys_L=["soft","medium","hard"]
    for CSC_Keys in CSC_Keys_L:
        Bands_Str_L_Reduced.append(CSC_Keys)
    print("Bands_Str_L_Reduced: ", Bands_Str_L_Reduced)
    #Key_L=["soft"]
    Key_L=Bands_Str_L_Reduced
    for Key in Key_L:
        Glob_Str="/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_*"+Key+".csv" #For Testing
        print("Glob_Str: ", Glob_Str)
        Path_L_Raw=glob.glob(Glob_Str)
        Path_L=[]
        for Path in Path_L_Raw:
            Src_Num=Find_Source_Num_Str(Path)
            if Src_Num=="All":
                continue
            else:
                Path_L.append(Path)
        Path_L.sort(key=Find_Source_Num)
        print("Path_L: ", Path_L)
        print("len(Path_L): ", len(Path_L))
        Data=pd.read_csv(Path_L[0])
        if(len(Path_L)>1):
            for i in range(1,len(Path_L)):
                Path=Path_L[i]
                Cur_Data=pd.read_csv(Path)
                Data=pd.concat([Data, Cur_Data])
        print("Data: ", Data)
        Outpath="/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_All_"+Key+".csv" #For Testing
        print("Outpath: ", Outpath)
        Data_CSV=Data.to_csv(Outpath,index=False)
        #print("Data_CSV: ", Data_CSV
        #COUNT_RATE
        #print("Data['COUNT_RATE']: ", Data["COUNT_RATE"])


#Flux_File_to_CSV_Big_Input(10125,9)
#Flux_File_to_CSV_Big_Input(10125,9,Key="soft")
#Flux_File_to_CSV_Big_Input(10125)
#Flux_File_to_CSV_Big_Input(10125,Key="soft")
#Flux_File_to_CSV_Big_Input("10125")
Merge_CSVs(10125)
#Flux_File_to_CSV_Big_Input_L([10125])
