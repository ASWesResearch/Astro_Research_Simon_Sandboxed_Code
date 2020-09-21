"""
This code is the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
"Exploring the High Energy Sky: A Survey of X-ray Sources in Galaxies Within 15 Megaparsecs of the Milky Way".
This thesis was submitted to the faculty of Wesleyan University in partial fulfillment of the requirements for the Degree of Bachelor of Arts
with Departmental Honors in Astronomy.

This version of the code is now modifed by Anthony Santini (Wesleyan 2019 Master's alumnus) for use in current research.


"""
import re
import os
import numpy as np
import glob
import subprocess as sp
##from done_email import when_done #Ant: This is no longer required
from ciao_contrib.runtool import *
from region import *
import csv

# Extracts counts from the from_csc sources processed by the raytracer and corrects for area. Creates counts_info.csv

def clean_results(obsid): #Ant: I don't know if this is still required for the refactored version of this code
    #the following makes all obsids into integers and removes anything globbed that isnt a number (and therefore not an obsid)
    for i in range(0,len(obsid)):
        try:
            obsid[i] = int(obsid[i])
        except ValueError:
            obsid[i] = -1
    obsid.remove(-1)
    obsid.remove(-1) #THIS IS A PRETTY POOR SOLUTION, BUT IT WORKS FOR NOW
    return obsid
'''
#Ant: This runs the specextract routine as a test. I don't know why it is still enabled but I am disabling it.
def spec(obsid):
    #Not used in this program -- can probably be deleted
    print "Beginning specextract for obsid " + str(obsid)
    obsid = str(obsid)
    print "Getting name of evt2 file..."
    ##name = glob.glob("chandra_from_csc/" + obsid + "/primary/*_evt2.fits")
    name = glob.glob("/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/*_evt2.fit*")
    print name
    print "Shortening name..."
    ##name = name[0][-24:] #Ant: What does this do? I'm waiting to see if I have to modify it
    print name
    print "Unlearning specextract..."
    os.system("punlearn specextract")
    print "Getting system path..."
    path = '/Volumes/xray/spirals/trace/' + obsid + '/'
    print "Getting number of files in the current obsid's directory..."
    num_files = len([f for f in os.listdir(path)
                if os.path.isfile(os.path.join(path, f))])
    print str(num_files) + " files found..."
    num = int((num_files-4)/5)
    print str(num) + " sources to be extracted..."
    for i in range(1,num+1):
        print "Running specextract..."
        print "Source " + str(i) + " of " + str(num) + " for obsid " + str(obsid) + "..."
        os.system("specextract \
        infile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + ".reg)]' \
        bkgfile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + "_bkg.reg)]' \
        outroot='chandra_from_csc/" + str(obsid) + "/primary/extracted/extracted_3747_" + str(i) + "' \
        correctpsf=yes \
        weight=no \
        clobber=yes \
        verbose=1")
'''
def effectiveAreaCorrect(count, cycle, hardness, instrument):
    """
    Corrects the input number of counts depending on the cycle and the wavelength range
    Count is an integer or float
    cycle is an integer 1-18
    hardness is 0 for .3-1 kev, 1 for 1-2.1 kev, 2 for 2.1-7.5 kev
    returns the corrected count
    Also need to somehow include acis s or i
    """

    print "Corecting count..."
    print "Original count: " + str(count)

    count = float(count)

    #######################################
    ###### Effective Area Correction ######
    #######################################

    s3EffCorrection = [1.0, 1.0, 1.0]
    s4EffCorrection = [1.17454, 1.01959, 1.01215]
    s5EffCorrection = [1.89781, 1.13248, 1.03366]
    s6EffCorrection = [1.75390, 1.11750, 1.02814]
    s7EffCorrection = [1.60477, 1.06976, 0.987487]
    s8EffCorrection = [1.53283, 1.01505, 0.926380]
    s9EffCorrection = [1.53859, 1.01607, 0.917762]
    s10EffCorrection = [1.54243, 1.01572, 0.915449]
    s11EffCorrection = [1.69040, 1.11287, 0.950149]
    s12EffCorrection = [1.81436, 1.07863, 0.948331]
    s13EffCorrection = [1.83137, 1.08009, 0.949277]
    s14EffCorrection = [1.84606, 1.08249, 0.952075]
    s15EffCorrection = [2.28854, 1.12162, 0.970426]
    s16EffCorrection = [2.79313, 1.21774, 0.978251]
    s17EffCorrection = [3.84065, 1.29621, 0.982907]
    s18EffCorrection = [4.41270, 1.33605, 0.990264]
    s19EffCorrection = [7.06738, 1.48238, 0.998766]

    i3EffCorrection = [1.0, 1.0, 1.0]
    i4EffCorrection = [1.00046, 1.01615, 1.01510]
    i5EffCorrection = [1.62312, 1.11492, 1.02181]
    i6EffCorrection = [1.51694, 1.10034, 1.01668]
    i7EffCorrection = [1.60943, 1.12987, 0.970840]
    i8EffCorrection = [1.52124, 1.06871, 0.909019]
    i9EffCorrection = [1.52345, 1.07181, 0.909649]
    i10EffCorrection = [1.58361, 1.11358, 0.943991]
    i11EffCorrection = [1.73488, 1.21974, 0.975594]
    i12EffCorrection = [1.75761, 1.19073, 0.975465]
    i13EffCorrection = [1.89543, 1.19287, 0.975131]
    i14EffCorrection = [1.91177, 1.19443, 0.975280]
    i15EffCorrection = [2.22792, 1.23182, 0.988722]
    i16EffCorrection = [2.68469, 1.33745, 0.995931]
    i17EffCorrection = [3.53210, 1.42364, 1.00066]
    i18EffCorrection = [3.99457, 1.46494, 1.00261]
    i19EffCorrection = [8.58820, 1.77972, 1.01788]

    sSoftAreaCorrection = [s3EffCorrection[0], s4EffCorrection[0], s5EffCorrection[0], s6EffCorrection[0], s7EffCorrection[0], s8EffCorrection[0], s9EffCorrection[0], s10EffCorrection[0], s11EffCorrection[0], s12EffCorrection[0], s13EffCorrection[0], s14EffCorrection[0], s15EffCorrection[0], s16EffCorrection[0], s17EffCorrection[0], s18EffCorrection[0], s19EffCorrection[0]]

    sMedAreaCorrection = [s3EffCorrection[1], s4EffCorrection[1], s5EffCorrection[1], s6EffCorrection[1], s7EffCorrection[1], s8EffCorrection[1], s9EffCorrection[1], s10EffCorrection[1], s11EffCorrection[1], s12EffCorrection[1], s13EffCorrection[1], s14EffCorrection[1], s15EffCorrection[1], s16EffCorrection[1], s17EffCorrection[1], s18EffCorrection[1], s19EffCorrection[1]]

    sHardAreaCorrection = [s3EffCorrection[2], s4EffCorrection[2], s5EffCorrection[2], s6EffCorrection[2], s7EffCorrection[2], s8EffCorrection[2], s9EffCorrection[2], s10EffCorrection[2], s11EffCorrection[2], s12EffCorrection[2], s13EffCorrection[2], s14EffCorrection[2], s15EffCorrection[2], s16EffCorrection[2], s17EffCorrection[2], s18EffCorrection[2], s19EffCorrection[2]]

    iSoftAreaCorrection = [i3EffCorrection[0], i4EffCorrection[0], i5EffCorrection[0], i6EffCorrection[0], i7EffCorrection[0], i8EffCorrection[0], i9EffCorrection[0], i10EffCorrection[0], i11EffCorrection[0], i12EffCorrection[0], i13EffCorrection[0], i14EffCorrection[0], i15EffCorrection[0], i16EffCorrection[0], i17EffCorrection[0], i18EffCorrection[0], i19EffCorrection[0]]

    iMedAreaCorrection = [i3EffCorrection[1], i4EffCorrection[1], i5EffCorrection[1], i6EffCorrection[1], i7EffCorrection[1], i8EffCorrection[1], i9EffCorrection[1], i10EffCorrection[1], i11EffCorrection[1], i12EffCorrection[1], i13EffCorrection[1], i14EffCorrection[1], i15EffCorrection[1], i16EffCorrection[1], i17EffCorrection[1], i18EffCorrection[1], i19EffCorrection[1]]

    iHardAreaCorrection = [i3EffCorrection[2], i4EffCorrection[2], i5EffCorrection[2], i6EffCorrection[2], i7EffCorrection[2], i8EffCorrection[2], i9EffCorrection[2], i10EffCorrection[2], i11EffCorrection[2], i12EffCorrection[2], i13EffCorrection[2], i14EffCorrection[2], i15EffCorrection[2], i16EffCorrection[2], i17EffCorrection[2], i18EffCorrection[2], i19EffCorrection[2]]

    val = cycle-3
    countOut = count

    if instrument == 's':
        if hardness == 0:
            if cycle == 1 or cycle == 2 or cycle == 3:
                countOut = count * sSoftAreaCorrection[0]
            else:
                countOut = count * sSoftAreaCorrection[val]
        if hardness == 1:
            if cycle == 1 or cycle == 2 or cycle == 3:
                countOut = count * sMedAreaCorrection[0]
            else:
                countOut = count * sMedAreaCorrection[val]
        if hardness == 2:
            if cycle == 1 or cycle == 2 or cycle == 3:
                countOut = count * sHardAreaCorrection[0]
            else:
                countOut = count * sHardAreaCorrection[val]
    elif instrument == 'i':
        if hardness == 0:
            if cycle == 1 or cycle == 2 or cycle == 3:
                countOut = count * iSoftAreaCorrection[0]
            else:
                countOut = count * iSoftAreaCorrection[val]
        if hardness == 1:
            if cycle == 1 or cycle == 2 or cycle == 3:
                countOut = count * iMedAreaCorrection[0]
            else:
                countOut = count * iMedAreaCorrection[val]
        if hardness == 2:
            if cycle == 1 or cycle == 2 or cycle == 3:
                countOut = count * iHardAreaCorrection[0]
            else:
                countOut = count * iHardAreaCorrection[val]

    print "Updated count: " + str(countOut)

    return countOut

def src_flux(obsid):
    """
    Takes as an input an obsid in the chandra_from_csc folder and extracts counts information for all sources
    """

    ###########################################################
    ###### Get the Number of Sources for the Given Obsid ######
    ###########################################################


    obsid = str(obsid)
    print "Beginning flux extraction for obsid " + obsid + "..."
    print "Getting name of evt2 file..."
    ##name = glob.glob("chandra_from_csc/" + obsid + "/primary/*_evt2.fits")
    Evt2_Fpath_L = glob.glob("/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/*_evt2.fit*")
    Evt2_Fpath=Evt2_Fpath_L[0]
    print "Evt2_Fpath: ", Evt2_Fpath
    ##print "Shortening name..."
    #Ant: This is done only for an exact directory on a computer not used anymore. It must be modified to create the shortened filename from an arbitary filepath
    ##name = name[0][-24:]
    name_L=Evt2_Fpath.split("/")
    print "name_L: ", name_L
    name=name_L[len(name_L)-1]
    print name
    print "Getting system path..."
    '''
    ##path = '/Volumes/xray/spirals/trace/' + obsid + '/' #Ant: This needs to be modfied to be using the Hybrid Regions files instead of the traced region files in the Trace directory
    #Ant: The purpose of this section is to count the total amount of souces for a given ObsID. This can be done better simply by counting the number of rows in the Hybrid_Region file for a given ObsID
    print "Getting number of files in the current obsid's directory..."
    num_files = len([f for f in os.listdir(path)
                if os.path.isfile(os.path.join(path, f))])
    print str(num_files) + " files found..."
    num = int((num_files-4)/5)
    print str(num) + " sources found..."
    '''
    #Ant: Caclulating the number of souces in the current observation
    Hybrid_Filepath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(obsid)+"/"+str(obsid)+"_Nearest_Neighbor_Hybrid.reg"
    Hybrid_File=open(Hybrid_Filepath)
    Hybrid_Region_Str=Hybrid_File.read()
    Header_String='# Region file format: DS9 version 3.0\nglobal color=blue font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0\n'
    Hybrid_Region_Str_L=Hybrid_Region_Str.split(Header_String)
    Hybrid_Region_Str_Reduced=Hybrid_Region_Str_L[1]
    #print "Hybrid_Region_Str_Reduced:\n", Hybrid_Region_Str_Reduced
    Hybrid_Region_Str_Reduced_L=Hybrid_Region_Str_Reduced.split("\n")
    #print "Hybrid_Region_Str_Reduced_L: ", Hybrid_Region_Str_Reduced_L
    #print "len(Hybrid_Region_Str_Reduced_L) Before Pop: ", len(Hybrid_Region_Str_Reduced_L)
    Hybrid_Region_Str_Reduced_L.pop(len(Hybrid_Region_Str_Reduced_L)-1)
    num=len(Hybrid_Region_Str_Reduced_L)
    print "num: ", num
    #######################################################
    ###### Get the Cycle and Instrument for the Data ######
    #######################################################

    print "Getting cycle..."
    ##date = sp.check_output("dmkeypar /Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + name + " DATE-OBS echo+", shell=True) #Ant: This need to be modifed to be used for all observations not just observations in CSC
    date = sp.check_output("dmkeypar /Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + name + " DATE-OBS echo+", shell=True)
    year = date[:4]
    shYear = year[2:]
    cycle = int(shYear)+1                           #eg. 2004 is cycle 5, 2005 is cycle 6, etc. This is pretty consistently true
    print "Got cycle " + str(cycle) + "..."


    # Will need to determine instrument below, not in this section, since it will change depending on which source we are looking at.

    #print "Determing instrument..."
    #currently do not know how to do this: assuming s for now to see if the rest of it works
    #in order to figure out what instrument the source is on, use the x,y coordinates of the source and get the location
    #inst = 's'


    ##################################
    ###### Initialize the Lists ######
    ##################################

    obs = []
    src = []
    rawCountsSoft = []
    rawCountsMed = []
    rawCountsHard = []
    areaSoft = []
    areaMed = []
    areaHard = []
    bkgCountsSoft = []
    bkgCountsMed = []
    bkgCountsHard = []
    bkgAreaSoft = []
    bkgAreaMed = []
    bkgAreaHard = []

    """
    Hybrid_BG_Filepath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(obsid)+"/"+str(obsid)+"_Nearest_Neighbor_Hybrid_Background.reg"
    Hybrid_BG_File=open(Hybrid_BG_Filepath)
    Hybrid_BG_Region_Str=Hybrid_BG_File.read()
    Header_String='# Region file format: DS9 version 3.0\nglobal color=blue font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0\n'
    Hybrid_BG_Region_Str_L=Hybrid_BG_Region_Str.split(Header_String)
    Hybrid_BG_Region_Str_Reduced=Hybrid_BG_Region_Str_L[1]
    #print "Hybrid_BG_Region_Str_Reduced:\n", Hybrid_BG_Region_Str_Reduced
    Hybrid_BG_Region_Str_Reduced_L=Hybrid_BG_Region_Str_Reduced.split("\n")
    #print "Hybrid_BG_Region_Str_Reduced_L: ", Hybrid_BG_Region_Str_Reduced_L
    #print "len(Hybrid_BG_Region_Str_Reduced_L) Before Pop: ", len(Hybrid_BG_Region_Str_Reduced_L)
    Hybrid_BG_Region_Str_Reduced_L.pop(len(Hybrid_BG_Region_Str_Reduced_L)-1)
    """
    #for Hybrid_BG_Region_Str in Hybrid_BG_Region_Str_Reduced_L: #Ant: I need to modify this so that it itterates through the Hybrid source background regions as well!
    for i in range(0,num):
        Hybrid_Region_Str=Hybrid_Region_Str_Reduced_L[i]
        print "Hybrid_Region_Str: ", Hybrid_Region_Str
        Hybrid_Region_Str_L=re.split("[(),]", Hybrid_Region_Str)
        print "Hybrid_Region_Str_L: ", Hybrid_Region_Str_L
        x_str=Hybrid_Region_Str_L[1]
        x=float(x_str)
        y_str=Hybrid_Region_Str_L[2]
        y=float(y_str)
        print "Got x and y values..."
        print "punlearn dmcoords..."
        os.system("punlearn dmcoords")
        print "running dmcoords..."
        #print "dmcoords '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "' option=sky x=" + str(x) + " y=" + str(y) + " celfmt=deg"
        os.system("dmcoords '/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "' option=sky x=" + str(x) + " y=" + str(y) + " celfmt=deg")
        print "Getting Chip ID..."
        chip = sp.check_output("pget dmcoords chip_id", shell=True)
        chip = int(chip)
        print "Got Chip ID: " + str(chip)

        if chip == 5 or chip == 7:
            inst = 's'
        else:
            inst = 'i'

        print inst

        #Ant: Reducing the Hybrid_Region_Str to a region shape
        print "Extracting soft X-rays for source " + str(i+1) + " in obsid " + obsid + "..."
        Hybrid_Region_Str_Reduced_Split_L=re.split("[; ]",Hybrid_Region_Str)
        print "Hybrid_Region_Str_Reduced_Split_L: ", Hybrid_Region_Str_Reduced_Split_L
        Hybrid_Region_Str_Reduced=Hybrid_Region_Str_Reduced_Split_L[1]
        print "Hybrid_Region_Str_Reduced: ", Hybrid_Region_Str_Reduced
        #Ant: Parsing the parsing the Hybrid source background region file to in order to feed that input into specextract
        """
        print "i: ", i
        j=i
        print "j Before: ", j
        """
        """
        if(j%2!=0):
            #j=j+1 #Ant: I lost a week to this bug...
            j=2*j
        """
        """
        j=2*j
        print "j After: ", j
        Hybrid_BG_Region_Outer_R_Str=Hybrid_BG_Region_Str_Reduced_L[j]
        Hybrid_BG_Region_Outer_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Outer_R_Str)
        Hybrid_BG_Region_Outer_R_Str_Reduced=Hybrid_BG_Region_Outer_R_Str_Reduced_L[1]
        print "Hybrid_BG_Region_Outer_R_Str_Reduced: ", Hybrid_BG_Region_Outer_R_Str_Reduced
        Hybrid_BG_Region_Inner_R_Str=Hybrid_BG_Region_Str_Reduced_L[j+1]
        Hybrid_BG_Region_Inner_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Inner_R_Str)
        Hybrid_BG_Region_Inner_R_Str_Reduced=Hybrid_BG_Region_Inner_R_Str_Reduced_L[1]
        print "Hybrid_BG_Region_Inner_R_Str_Reduced: ", Hybrid_BG_Region_Inner_R_Str_Reduced
        Hybrid_BG_Region=Hybrid_BG_Region_Outer_R_Str_Reduced+Hybrid_BG_Region_Inner_R_Str_Reduced
        """
        print "i: ", i
        Source_Num=i+1
        print "Source_Num: ", Source_Num
        #/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/10125/Individual_Source_Regions/Source_1_ObsID_10125_Nearest_Neighbor_Hybrid_Overlap_Corrected_Background.reg
        Hybrid_BG_Filepath="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"+str(obsid)+"/Individual_Source_Regions/Source_"+str(Source_Num)+"_ObsID_"+str(obsid)+"_Nearest_Neighbor_Hybrid_Overlap_Corrected_Background.reg"
        print "Hybrid_BG_Filepath: ", Hybrid_BG_Filepath
        #print "Hybrid_BG_Region: ", Hybrid_BG_Region

        Hybrid_BG_File=open(Hybrid_BG_Filepath)
        Hybrid_BG_Region_Str=Hybrid_BG_File.read()
        Hybrid_BG_File.close()
        Header_String='# Region file format: DS9 version 3.0\nglobal color=blue font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0\n'
        Hybrid_BG_Region_Str_L=Hybrid_BG_Region_Str.split(Header_String)
        Hybrid_BG_Region_Str_Reduced=Hybrid_BG_Region_Str_L[1]
        #print "Hybrid_BG_Region_Str_Reduced:\n", Hybrid_BG_Region_Str_Reduced
        Hybrid_BG_Region_Str_Reduced_L=Hybrid_BG_Region_Str_Reduced.split("\n")
        #print "Hybrid_BG_Region_Str_Reduced_L: ", Hybrid_BG_Region_Str_Reduced_L
        #print "len(Hybrid_BG_Region_Str_Reduced_L) Before Pop: ", len(Hybrid_BG_Region_Str_Reduced_L)
        Hybrid_BG_Region_Str_Reduced_L.pop(len(Hybrid_BG_Region_Str_Reduced_L)-1)
        Hybrid_BG_Region=""
        for Region in Hybrid_BG_Region_Str_Reduced_L:
            Region_Reduced=re.split("[; ]",Region)[1]
            Hybrid_BG_Region=Hybrid_BG_Region+Region_Reduced
        #Ant: Calculating the area of the current Hybrid source region
        #shape1 ='circle(' + str(X_Phys) +','+ str(Y_Phys)+','+ str(cur_r)+')' #shape1:-str, shape1, The shape string of the current area circle
        shape1=Hybrid_Region_Str_Reduced
        #print "shape1 : ",shape1
        r1 = regParse(shape1) #r1:-Region, Region 1, the region of the current area circle
        Cur_Source_Area = regArea(r1,0,0,8192,8192,1) #a1_cur:-float, Area_1_Current, The area of the current area circle
        shape2=Hybrid_BG_Region
        r2= regParse(shape2)
        Cur_Background_Area= regArea(r2,0,0,8192,8192,1)
        directory="/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/"
        if not os.path.exists(directory):
            os.makedirs(directory)

        #########################
        ###### Soft X-rays ######
        #########################

        """
        os.system("dmextract \
        infile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=300:1000][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + ".reg)]\" \
        outfile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_softcounts.fits\" \ #Ant: This outfile needs to be modifed to be in the Simon Sandboxed Code Directory
        opt=generic \
        bkg=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=300:1000][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + "_bkg.reg)]\" \
        clobber=no")
        """
        os.system("dmextract \
        infile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[energy=300:1000][sky="+Hybrid_Region_Str_Reduced+"]\" \
        outfile=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_softcounts.fits\" \
        opt=generic \
        bkg=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[energy=300:1000][sky=region("+str(Hybrid_BG_Filepath)+")]\" \
        clobber=no")

        #Getting counts
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_softcounts.fits[cols counts]'")
        count = sp.check_output("pget dmstat out_sum", shell=True)
        count = effectiveAreaCorrect(count, cycle, 0, inst)
        rawCountsSoft.append(count)

        """
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_softcounts.fits[cols area]'") #Ant: This isn't working. It needs to be replaced with ciao region area to calculate geometeric area
        area = sp.check_output("pget dmstat out_sum", shell=True)
        areaSoft.append(area)
        """
        areaSoft.append(Cur_Source_Area)

        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_softcounts.fits[cols bg_counts]'")
        bgCount = sp.check_output("pget dmstat out_sum", shell=True)
        bgCount = effectiveAreaCorrect(bgCount, cycle, 0, inst)
        bkgCountsSoft.append(bgCount)

        """
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_softcounts.fits[cols bg_area]'")
        bgArea = sp.check_output("pget dmstat out_sum", shell=True)
        bkgAreaSoft.append(bgArea)
        """
        bkgAreaSoft.append(Cur_Background_Area)

        ########################
        ###### Med X-rays ######
        ########################

        print "Extracting medium X-rays for source " + str(i+1) + " in obsid " + obsid + "..."

        os.system("dmextract \
        infile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[energy=1000:2100][sky="+Hybrid_Region_Str_Reduced+"]\" \
        outfile=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_medcounts.fits\" \
        opt=generic \
        bkg=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[energy=1000:2100][sky=region("+str(Hybrid_BG_Filepath)+")]\" \
        clobber=no")

        #Getting counts
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_medcounts.fits[cols counts]'")
        count = sp.check_output("pget dmstat out_sum", shell=True)
        count = effectiveAreaCorrect(count, cycle, 1, inst)
        rawCountsMed.append(count)

        """
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_medcounts.fits[cols area]'") #Ant: This isn't working. It needs to be replaced with ciao region area to calculate geometeric area
        area = sp.check_output("pget dmstat out_sum", shell=True)
        areaMed.append(area)
        """
        areaMed.append(Cur_Source_Area)

        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_medcounts.fits[cols bg_counts]'")
        bgCount = sp.check_output("pget dmstat out_sum", shell=True)
        bgCount = effectiveAreaCorrect(bgCount, cycle, 1, inst)
        bkgCountsMed.append(bgCount)

        """
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_medcounts.fits[cols bg_area]'")
        bgArea = sp.check_output("pget dmstat out_sum", shell=True)
        bkgAreaMed.append(bgArea)
        """
        bkgAreaMed.append(Cur_Background_Area)

        #########################
        ###### Hard X-rays ######
        #########################

        print "Extracting hard X-rays for source " + str(i+1) + " in obsid " + obsid + "..."
        #Ant: This needs to be updated to use the Hybrid region files
        os.system("dmextract \
        infile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[energy=2100:7500][sky="+Hybrid_Region_Str_Reduced+"]\" \
        outfile=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_hardcounts.fits\" \
        opt=generic \
        bkg=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[energy=2100:7500][sky=region("+str(Hybrid_BG_Filepath)+")]\" \
        clobber=no")

        #Getting counts
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_hardcounts.fits[cols counts]'")
        count = sp.check_output("pget dmstat out_sum", shell=True)
        count = effectiveAreaCorrect(count, cycle, 2, inst)
        rawCountsHard.append(count)

        """
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_hardcounts.fits[cols area]'") #Ant: This isn't working. It needs to be replaced with ciao region area to calculate geometeric area
        area = sp.check_output("pget dmstat out_sum", shell=True)
        areaHard.append(area)
        """
        areaHard.append(Cur_Source_Area)

        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_hardcounts.fits[cols bg_counts]'")
        bgCount = sp.check_output("pget dmstat out_sum", shell=True)
        bgCount = effectiveAreaCorrect(bgCount, cycle, 2, inst)
        bkgCountsHard.append(bgCount)

        """
        os.system("dmstat '/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/extracted_counts_info/"+obsid+"/" + obsid + "_" + str(i+1) + "_hardcounts.fits[cols bg_area]'")
        bgArea = sp.check_output("pget dmstat out_sum", shell=True)
        bkgAreaHard.append(bgArea)
        """
        bkgAreaHard.append(Cur_Background_Area)

        obs.append(obsid)
        src.append(i+1)

        #Reset the values - unsure if this is necessary
        # count = 0
        # area = 0
        # bgCount = 0
        # bgArea = 0

    # print len(obs)
    # print len(src)
    # print len(rawCountsSoft)
    # print len(rawCountsMed)
    # print len(rawCountsHard)
    # print len(areaSoft)
    # print len(areaMed)
    # print len(areaHard)
    # print len(bkgCountsSoft)
    # print len(bkgCountsMed)
    # print len(bkgCountsHard)
    # print len(bkgAreaSoft)
    # print len(bkgAreaMed)
    # print len(bkgAreaHard)

    #############################
    ###### Return the Data ######
    #############################

    return obs, src, rawCountsSoft, rawCountsMed, rawCountsHard, areaSoft, areaMed, areaHard, bkgCountsSoft, bkgCountsMed, bkgCountsHard, bkgAreaSoft, bkgAreaMed, bkgAreaHard



##if __name__ == '__main__': #Ant: I don't know if this is needed but I am quoting it out anyway
def Driver(obsid_L):
    """
    print "Reading in obsids... "
    obsid_path = glob.glob("chandra_from_csc/*") #Ant: This should be modifed to just use an input list of ObsIds

    obsid = []
    for i in range(0,len(obsid_path)):
        #print obsid_path[i]
        obsid.append(obsid_path[i][17:])
    print "...done!"


    obsid = clean_results(obsid) #Ant: I don't know if this is still required for the refactored version of this code
    """

    ##############################
    ###### Start Processing ######
    ##############################


    #print "All obsid names processed."

    obsidData = ["OBSID"]
    sourceData = ["SOURCE"]
    rawCountsSoft = ["RAW_COUNTS[.3-1]"]
    rawCountsMed = ["RAW_COUNTS[1-2.1]"]
    rawCountsHard = ["RAW_COUNTS[2.1-7.5]"]
    #areaSoft = ["AREA[.3-1]"]
    areaSoft = ["AREA"]
    areaMed = ["AREA[1-2.1]"]
    areaHard = ["AREA[2.1-7.5]"]
    bkgCountsSoft = ["BKG_COUNTS[.3-1]"]
    bkgCountsMed = ["BKG_COUNTS[1-2.1]"]
    bkgCountsHard = ["BKG_COUNTS[2.1-7.5]"]
    bkgAreaSoft = ["BKG_AREA"]
    bkgAreaMed = ["BKG_AREA[1-2.1]"]
    bkgAreaHard = ["BKG_AREA[2.1-7.5]"]


    for i in range(0,len(obsid_L)):
        try: #Ant: ALL ERRORS SHOULD NOT BE ASSUMED TO BE THE RESULT OF A FAILED FLUX EXTRACTION! THIS NEEDS TO BE FIXED AND TAKEN OUT OF THE TRY/EXCEPT CONDITIONAL ENTIRELY!
            print "\nAttempting flux retrieval..."
            obs, src, countS, countM, countH, areaS, areaM, areaH, bgcountS, bgcountM, bgcountH, bgareaS, bgareaM, bgareaH = src_flux(str(obsid_L[i]))

            # print obs
            # print src
            # print countS
            # print countM
            # print countH
            # print areaS
            # print areaM
            # print areaH
            # print bgcountS
            # print bgcountM
            # print bgcountH
            # print bgareaS
            # print bgareaM
            # print bgareaH


            for j in range(0,len(obs)):
                obsidData.append(obs[j])
                sourceData.append(src[j])
                rawCountsSoft.append(countS[j])
                rawCountsMed.append(countM[j])
                rawCountsHard.append(countH[j])
                #areaSoft.append(areaS[j][:-1])
                areaSoft.append(areaS[j])
                #areaMed.append(areaM[j][:-1])
                areaMed.append(areaM[j])
                #areaHard.append(areaH[j][:-1])
                areaHard.append(areaH[j])
                bkgCountsSoft.append(bgcountS[j])
                bkgCountsMed.append(bgcountM[j])
                bkgCountsHard.append(bgcountH[j])
                #bkgAreaSoft.append(bgareaS[j][:-1])
                bkgAreaSoft.append(bgareaS[j])
                #bkgAreaMed.append(bgareaM[j][:-1])
                bkgAreaMed.append(bgareaM[j])
                #bkgAreaHard.append(bgareaH[j][:-1])
                bkgAreaHard.append(bgareaH[j])


        except OSError:
            #print "No traced sources in obsid " + str(obsid[i])
            print "\nFlux extraction failed for obsid " + str(obsid[i])

    #zippedList = zip(obsidData, sourceData, rawCountsSoft, rawCountsMed, rawCountsHard, areaSoft, areaMed, areaHard, bkgCountsSoft, bkgCountsMed, bkgCountsHard, bkgAreaSoft, bkgAreaMed, bkgAreaHard)
    zippedList = zip(obsidData, sourceData, rawCountsSoft, rawCountsMed, rawCountsHard, bkgCountsSoft, bkgCountsMed, bkgCountsHard, areaSoft, bkgAreaSoft)


    # print len(obsidData)
#
    # print obsidData
    # for i in range(0,len(obsidData)):
        # print zippedList[i]

    #with open('counts_info.csv', 'wb') as csvfile:
    with open('/Volumes/xray/anthony/Simon_Sandboxed_Code/Hardness_Ratios/counts_info.csv', 'wb') as csvfile:
        cWriter = csv.writer(csvfile, delimiter=',')
        for i in range(0,len(obsidData)):
            cWriter.writerow(zippedList[i])


    #when_done("specextract")

def Main():
    #Driver([10125])
    Driver([6096, 1971, 1972, 768, 952, 11674, 13255, 13253, 13246, 12952, 12953, 13247, 12951, 2025, 9548, 2149, 2197, 9510, 6131, 5908, 803, 14342, 12995, 2064, 16024, 12992, 14332, 13202, 793, 2933, 11104, 379, 2056, 2055, 2922, 9506, 11344, 766, 4688, 6869, 6872, 3554, 2057, 2058, 8041, 9121, 9546, 7252, 7060, 9553, 5930, 5931, 5929, 2079, 5905, 9527, 4689, 3947, 1563, 9507, 4613, 794, 11775, 11271, 3951, 2062, 2027, 2060, 2061, 2070, 2032, 7154, 7153, 11779, 5932, 2976, 4613, 794, 1043, 4632, 4631, 4633, 4404, 2059, 12095, 2040, 2915, 4372, 2069, 11229, 7848, 15383, 10125, 2031, 10875, 12889, 12888, 321, 322, 9551, 9550, 3954, 2020, 2068, 4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552])

Main()
