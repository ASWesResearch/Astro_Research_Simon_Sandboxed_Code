import os
import sys
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import glob

#Constants:
#Root_Path="/Volumes/"
Root_Path="/opt/"

print(sys.version)
sys.path.append(Root_Path+"xray/anthony/Research_Git/")
from ObsID_From_CSV_Query import ObsID_From_CSV_Query
#Exception Classes
class FileNotFoundException(Exception):
    "Raised when a file is not found"
    def __init__(self, ObsID):
        self.ObsID = ObsID
        self.message = "ObsID "+str(self.ObsID)+" is Missing File"
        super().__init__(self.message)

def Find_Source_Num(Fpath):
    Fname=Fpath.split("/")[-1]
    Src_Num_Str=Fname.split("_")[3]
    Src_Num=int(Src_Num_Str)
    return Src_Num

def Find_ObsID(Fpath):
    Fname=Fpath.split("/")[-1]
    ObsID_Str=Fname.split("_")[2]
    ObsID=int(ObsID_Str)
    return ObsID

def Find_Band(Fpath):
    Fname=Fpath.split("/")[-1]
    Band_Str_Raw=Fname.split("_")[-1]
    Band=Band_Str_Raw.split(".flux")[0]
    return Band

def Flux_File_to_CSV(Path, Outpath):
    with fits.open(Path, memmap=True) as hdu_list: # set memmap=True for large files
        hdu = hdu_list[1]
        table = Table(hdu.data)
        names = [name for name in table.colnames if len(table[name].shape) <= 1]
        #df = table.to_pandas()
        #tbl[names].to_pandas(...)
        df=table[names].to_pandas()
        #print("Data_Frame:\n", df)
        Band=Find_Band(Path)
        df.insert(loc = 0, column = 'Band', value = Band)
        ObsID=Find_ObsID(Path)
        df.insert(loc = 0, column = 'ObsID', value = ObsID)
        Src_Num=Find_Source_Num(Path)
        df.insert(loc = 0, column = 'Source_Num', value = Src_Num)
        print("Data_Frame:\n", df)
        #print(table[0])
        #table=table[0]
        #table.write('Test_Table_3.ecsv', delimiter='\t', format='ascii')
        #table.write('Test_Table_3.ecsv', delimiter=' ', format='ascii')
        df.to_csv(Outpath, index=False)

#Flux_File_to_CSV("/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/10125/Source_Flux_10125_1_soft.flux")

def Flux_File_to_CSV_Big_Input(ObsID="All",Source_Num="All",Key="All"):
    if(Key=="All"):
        Key="*"
    if(Source_Num=="All"):
        Source_Num="*"
    if(ObsID=="All"):
        ObsID="*"
    Glob_Str="/opt/xray/anthony/expansion_backup/Source_Fluxes/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_"+str(Source_Num)+"_"+Key+".flux" #Real One
    #/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/10125
    #Glob_Str="/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_"+str(Source_Num)+"_"+Key+".flux" #For Testing
    #/opt/xray/anthony/expansion_backup/Source_Fluxes/
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

def Find_Source_Num_Str(Fpath):
    Fname=Fpath.split("/")[-1]
    Src_Num_Str=Fname.split("_")[3]
    return Src_Num_Str

def Merge_CSVs(ObsID):
    Glob_Str="/opt/xray/anthony/expansion_backup/Source_Fluxes/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_*.csv"
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
        #Glob_Str="/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_*"+Key+".csv" #For Testing
        #/opt/xray/anthony/expansion_backup/Source_Fluxes/
        Glob_Str="/opt/xray/anthony/expansion_backup/Source_Fluxes/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_*"+Key+".csv" #Real One
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
        try:
            if(len(Path_L)==0):
                print("ObsID "+str(ObsID)+" is Missing File")
                raise FileNotFoundException(ObsID)
        except:
            print("ObsID "+str(ObsID)+" is Missing File")
            raise FileNotFoundException(ObsID)
        Data=pd.read_csv(Path_L[0])
        if(len(Path_L)>1):
            for i in range(1,len(Path_L)):
                Path=Path_L[i]
                Cur_Data=pd.read_csv(Path)
                Data=pd.concat([Data, Cur_Data])
        print("Data: ", Data)
        #Outpath="/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_All_"+Key+".csv" #For Testing
        #/opt/xray/anthony/expansion_backup/Source_Fluxes
        Outpath="/opt/xray/anthony/expansion_backup/Source_Fluxes/"+str(ObsID)+"/Source_Flux_"+str(ObsID)+"_All_"+Key+".csv"  #Real one
        print("Outpath: ", Outpath)
        Data_CSV=Data.to_csv(Outpath,index=False)
        #print("Data_CSV: ", Data_CSV
        #COUNT_RATE
        #print("Data['COUNT_RATE']: ", Data["COUNT_RATE"])

def Flux_File_to_CSV_Big_Input_L(ObsID_L, Clobber_Bool=False):
    Error_List=[]
    for ObsID in ObsID_L:
        if(Clobber_Bool==False):
            #/opt/xray/anthony/expansion_backup/Source_Fluxes/10125/*.csv
            Glob_Path="/opt/xray/anthony/expansion_backup/Source_Fluxes/"+str(ObsID)+"/*.csv"
            Glob_L=glob.glob(Glob_Path)
            if(len(Glob_L)>0):
                print(str(ObsID)+" Already exists and Clobber is set to False")
                continue
        try:
            Flux_File_to_CSV_Big_Input(ObsID=ObsID)
            Merge_CSVs(ObsID)
        except:
            Error_List.append(ObsID)
    print("Error_List: ", Error_List)
#Flux_File_to_CSV_Big_Input(10125,9)
#Flux_File_to_CSV_Big_Input(10125,9,Key="soft")
#Flux_File_to_CSV_Big_Input(10125)
#Flux_File_to_CSV_Big_Input(10125,Key="soft")
#Flux_File_to_CSV_Big_Input("10125")
#Merge_CSVs(10125)
#Flux_File_to_CSV_Big_Input_L([316, 318, 349, 350, 353, 354, 361, 378, 379, 383, 384, 388, 389, 390, 391, 392, 393, 394, 395, 402, 404, 405, 407, 409, 410, 411, 413, 414, 735, 782, 784, 790, 792, 793, 794, 795, 797, 808, 864, 870, 871, 872, 882, 934, 942, 962, 969, 1302, 1564, 1579, 1586, 1587, 1618, 1622, 1633, 1634, 1635, 1636, 1637, 1638, 1640, 2014, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2039, 2040, 2057, 2058, 2059, 2064, 2065, 2075, 2076, 2148, 2197, 2198, 2255, 2260, 2340, 2686, 2879, 2885, 2916, 2917, 2918, 2919, 2925, 2933, 2934, 2949, 2950, 2976, 2978, 3012, 3325, 3550, 3551, 3786, 3787, 3788, 3925, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 3937, 3938, 3939, 3940, 3941, 3942, 3943, 3949, 3950, 3953, 3954, 3965, 4010, 4016, 4017, 4019, 4176, 4555, 4556, 4557, 4558, 4613, 4627, 4628, 4629, 4630, 4688, 4689, 4690, 4692, 4693, 4694, 4696, 4697, 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 4734, 4735, 4736, 4737, 4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 5197, 5283, 5296, 5297, 5300, 5301, 5302, 5309, 5322, 5323, 5337, 5338, 5339, 5340, 5619, 5644, 5905, 5911, 5929, 5930, 5931, 5935, 5936, 5937, 5938, 5939, 5940, 5941, 5942, 5943, 5944, 5945, 5946, 5947, 5948, 5949, 6096, 6097, 6114, 6115, 6118, 6152, 6169, 6170, 6175, 6184, 6185, 6361, 6727, 6781, 6782, 7060, 7069, 7073, 7074, 7075, 7076, 7082, 7083, 7084, 7086, 7087, 7090, 7091, 7093, 7095, 7096, 7098, 7101, 7103, 7104, 7105, 7106, 7111, 7113, 7115, 7116, 7118, 7120, 7121, 7123, 7124, 7127, 7132, 7134, 7146, 7147, 7150, 7152, 7153, 7154, 7252, 7369, 7797, 7798, 7799, 7800, 7850, 7858, 7863, 7885, 8050, 8052, 8053, 8058, 8086, 8091, 8098, 8125, 8126, 8190, 8197, 8198, 8458, 8464, 8465, 8489, 8490, 9120, 9121, 9122, 9278, 9532, 9533, 9534, 9535, 9536, 9537, 9538, 9539, 9540, 9541, 9542, 9543, 9545, 9546, 9547, 9548, 9549, 9550, 9551, 9552, 9553, 9570, 9805, 9877, 9883, 10025, 10026, 10027, 10125, 10274, 10275, 10276, 10277, 10278, 10279, 10280, 10281, 10289, 10290, 10291, 10292, 10293, 10542, 10543, 10544, 10545, 10559, 10560, 10722, 10723, 10724, 10725, 10726, 10868, 10875, 10925, 11032, 11033, 11034, 11080, 11081, 11082, 11083, 11084, 11085, 11086, 11104, 11260, 11268, 11272, 11273, 11289, 11295, 11309, 11311, 11317, 11761, 11782, 11786, 11800, 11846, 11847, 12095, 12155, 12156, 12238, 12239, 12301, 12437, 12473, 12562, 12668, 12696, 12748, 12978, 12981, 12992, 12993, 12994, 12995, 12996, 13018, 13202, 13241, 13248, 13303, 13304, 13439, 13686, 13728, 13765, 13791, 13796, 13812, 13813, 13814, 13815, 13816, 13817, 13819, 13820, 13821, 13822, 13829, 13830, 13831, 13832, 14017, 14018, 14230, 14231, 14332, 14341, 14342, 14349, 14350, 14351, 14376, 14378, 14383, 14384, 14412, 14419, 14437, 14442, 14471, 14675, 14676, 14795, 14801, 14896, 14902, 14912, 14984, 14985, 15149, 15190, 15200, 15294, 15295, 15333, 15382, 15384, 15386, 15387, 15496, 15553, 15572, 15574, 15579, 15582, 15587, 15588, 15589, 15594, 15603, 15616, 15646, 15756, 15760, 15771, 15803, 16000, 16001, 16002, 16003, 16005, 16023, 16024, 16028, 16029, 16032, 16033, 16068, 16069, 16121, 16122, 16234, 16260, 16261, 16262, 16276, 16277, 16484, 16485, 16556, 16580, 16745, 16969, 16978, 16983, 16991, 16994, 16995, 16996, 16997, 17000, 17003, 17007, 17032, 17155, 17180, 17461, 17462, 17471, 17472, 17547, 17569, 17570, 17571, 17578, 17678, 17890, 17891, 18047, 18048, 18053, 18054, 18062, 18063, 18064, 18065, 18066, 18067, 18068, 18069, 18070, 18071, 18072, 18073, 18340, 18341, 18342, 18343, 18352, 18440, 18454, 18455, 18461, 18462, 18760, 18875, 19297, 19304, 19339, 19344, 19345, 19346, 19348, 19350, 19351, 19354, 19357, 19363, 19374, 19386, 19387, 19392, 19393, 19394, 19397, 19403, 19407, 19411, 19414, 19416, 19417, 19421, 19422, 19428, 19437, 19497, 19521, 19522, 19524, 19747, 19748, 19981, 19982, 20333, 20343, 20353, 20356, 20495, 20585, 20752, 20753, 20794, 20965, 20966, 20992, 20993, 20997, 20998, 20999, 21000, 21001, 21003, 21036, 21077, 21082, 21230, 21350, 21351, 21352, 21384, 21471, 21472, 21473, 21474, 21479, 21545, 21639, 21640, 21647, 21648, 21649, 21698, 21699, 21853, 22189, 22194, 22372, 22375, 22478, 22479, 22480, 22481, 22482, 22714, 22715, 23075, 23076, 23140, 23141, 23216, 23217, 23218, 23219, 23220, 23223, 23266, 23472, 23473, 23474, 23475, 23476, 23477, 23478, 23479, 23480, 23481, 23482, 23483, 23484, 23485, 23486, 23487, 23488, 23489, 23490, 23491, 23492, 23493, 23494, 23495, 23496, 23497, 23498, 23499, 23559, 23561, 23564, 24707, 24981, 24986, 25179, 25186, 25191, 25689]) #Full List with Unarchived ObsIDs 23500 and 23501 removed #Removed ObsIDs 380,400,963,1578 due to missing files. #Removing 1621, 1624 until issues with it are resolved.
#ObsID_L=ObsID_From_CSV_Query.Read_ObsIDs(Remove_Unarchived=True)
ObsID_L=ObsID_From_CSV_Query.Read_ObsIDs(Remove_Unarchived=True, ObsID_Tester_Bool=True)
Flux_File_to_CSV_Big_Input_L(ObsID_L)
