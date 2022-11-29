import pandas as pd
from astropy.io import fits
from astropy.table import Table

def Flux_File_to_CSV(Path):
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
        df.to_csv("Test_Table_3.csv")

Flux_File_to_CSV("/opt/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/10125/Source_Flux_10125_1_soft.flux")
