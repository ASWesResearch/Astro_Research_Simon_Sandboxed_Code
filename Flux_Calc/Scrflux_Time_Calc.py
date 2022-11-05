import matplotlib.pyplot as plt
import pandas as pd
def Srcflux_Plot(Path):
    Data=pd.read_csv(Path)
    #RAW_COUNTS[.3-1]	RAW_COUNTS[1-2.1]	RAW_COUNTS[2.1-7.5]
    Data["RAW_COUNTS[.3-7.5]"]=Data["RAW_COUNTS[.3-1]"]+Data["RAW_COUNTS[1-2.1]"]+Data["RAW_COUNTS[2.1-7.5]"]
    Data=Data.sort_values(by=["RAW_COUNTS[.3-7.5]"])
    print Data
    plt.semilogx(Data["RAW_COUNTS[.3-7.5]"],Data["Time"],".")
    plt.xlabel("Counts")
    plt.ylabel("Time(s)")
    #plt.show()
    plt.savefig("Time_vs_Counts.pdf")
    """
    plt.clf()
    Data["Hardness_Ratio"]=(Data["RAW_COUNTS[2.1-7.5]"]-Data["RAW_COUNTS[.3-1]"])/(Data["RAW_COUNTS[2.1-7.5]"]+Data["RAW_COUNTS[.3-1]"])
    Data=Data.sort_values(by=["Hardness_Ratio"])
    plt.plot(Data["Hardness_Ratio"],Data["Time"],".")
    #plt.show()
    #plt.savefig("Hardness_Ratio_vs_Counts.pdf")
    """
Srcflux_Plot("/Volumes/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Srcflux_Times.csv")
