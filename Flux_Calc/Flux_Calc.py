"""
This code is based on the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
"Exploring the High Energy Sky: A Survey of X-ray Sources in Galaxies Within 15 Megaparsecs of the Milky Way".
This thesis was submitted to the faculty of Wesleyan University in partial fulfillment of the requirements for the Degree of Bachelor of Arts
with Departmental Honors in Astronomy.

This version of the code is now modifed by Anthony Santini (Wesleyan 2019 Master's alumnus) for use in current research.


"""

import os
#from ciao_contrib.runtool import *
import numpy as np
import glob
import re

def spec(obsid):
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
    print "Unlearning specextract..."
    os.system("punlearn specextract")
    #print "Getting system path..."
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
    #for Hybrid_BG_Region_Str in Hybrid_BG_Region_Str_Reduced_L: #Ant: I need to modify this so that it itterates through the Hybrid source background regions as well!
    #for i in range(1,num+1):
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
        RA = sp.check_output("pget dmcoords ra", shell=True)
        Dec = sp.check_output("pget dmcoords dec", shell=True)
        print inst
        #Ant: Reducing the Hybrid_Region_Str to a region shape
        print "Extracting soft X-rays for source " + str(i+1) + " in obsid " + obsid + "..."
        Hybrid_Region_Str_Reduced_Split_L=re.split("[; ]",Hybrid_Region_Str)
        print "Hybrid_Region_Str_Reduced_Split_L: ", Hybrid_Region_Str_Reduced_Split_L
        Hybrid_Region_Str_Reduced=Hybrid_Region_Str_Reduced_Split_L[1]
        print "Hybrid_Region_Str_Reduced: ", Hybrid_Region_Str_Reduced
        #Ant: Parsing the parsing the Hybrid source background region file to in order to feed that input into specextract
        j=i
        if(j%2!=0):
            j=j+1
        Hybrid_BG_Region_Outer_R_Str=Hybrid_BG_Region_Str_Reduced_L[j]
        Hybrid_BG_Region_Outer_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Outer_R_Str)
        Hybrid_BG_Region_Outer_R_Str_Reduced=Hybrid_BG_Region_Outer_R_Str_Reduced_L[1]
        print "Hybrid_BG_Region_Outer_R_Str_Reduced: ", Hybrid_BG_Region_Outer_R_Str_Reduced
        Hybrid_BG_Region_Inner_R_Str=Hybrid_BG_Region_Str_Reduced_L[j+1]
        Hybrid_BG_Region_Inner_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Inner_R_Str)
        Hybrid_BG_Region_Inner_R_Str_Reduced=Hybrid_BG_Region_Inner_R_Str_Reduced_L[1]
        print "Hybrid_BG_Region_Inner_R_Str_Reduced: ", Hybrid_BG_Region_Inner_R_Str_Reduced
        Hybrid_BG_Region=Hybrid_BG_Region_Outer_R_Str_Reduced+Hybrid_BG_Region_Inner_R_Str_Reduced

        # print "specextract " + obsid + " " + str(i)
        # print "Assigning new defaults to specextract for source number " + str(i) + "..."
        # print "Assigning infile..."
        # os.system("pset specextract infile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + ".reg)]'")
        # print "Assigning bkgfile..."
        # os.system("pset specextract bkgfile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + "_bkg.reg)]'")
        # print "Assigning output file..."
        # os.system("pset specextract outroot='chandra_from_csc/" + str(obsid) + "/primary/extracted/extracted_3747_" + str(i) + "'")
        # print "Correct psf..."
        # os.system("pset specextract correctpsf=yes")
        # print "Weight..."
        # os.system("pset specextract weight=no")
        # print "Running specextract..."
        # os.system("specextract clobber=yes verbose=5")
        # print "...done with source number " + str(i) + "!"

        #ALTERNATELY (USE THIS INSTEAD SO SPECEXTRACT DOES NOT BUG YOU FOR CONFIRMATIONS):
        print "Running specextract..."
        print "Source " + str(i+1) + " of " + str(num) + " for obsid " + str(obsid) + "..."
        """
        os.system("specextract \
        infile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + ".reg)]' \
        bkgfile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + "_bkg.reg)]' \
        outroot='chandra_from_csc/" + str(obsid) + "/primary/extracted/extracted_" + obsid + "_" + str(i) + "' \
        correctpsf=yes \
        weight=no \
        clobber=yes \
        verbose=1")
        """
        #Ant: The string syntax for this may be incorrect

        os.system("specextract \
        infile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_Region_Str_Reduced+"]\" \
        bkgfile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_BG_Region+"]\" \
        outroot=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/" + str(obsid) +"/"+ "extracted_spectra" + obsid + "_" + str(i+1) + "' \
        correctpsf=yes \
        weight=no \
        clobber=yes \
        verbose=1")

        #Ant: Need to check string syntax here!
        #Ant: Need pass the input files from Specextract! I don't know what exact ones to pass. Should I pass all of them?
        os.system("srcflux \
        infile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) +"\" \
        pos=\"RA+" "+Dec\
        outroot=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes/" + str(obsid) +"/"+ "Source_Flux_" + obsid + "_" + str(i+1) + "\
        bands=\"0.3:7.0:1.9 " "\
        srcreg=Hybrid_Region_Str_Reduced\
        bkgreg=Hybrid_BG_Region \
        psfmethod=ideal \
        arffile=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".arf" "\
        rmffile=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".rmf" "\
        clobber=no \
        verbose=1")

##if __name__ == '__main__':
def Driver(obsid_L):
    for i in range(0,len(obsid_L)):
        try:
            #print "Attempting specextract..."
            spec(obsid_L[i])
        except OSError:
            print "No traced sources in obsid " + str(obsid[i])
            print "specextract failed for obsid " + str(obsid[i])

    #when_done("specextract") #Ant: Disabled to not accidently email Simon

def Main():
    #Driver([10125])
    Driver([6096, 1971, 1972, 768, 952, 11674, 13255, 13253, 13246, 12952, 12953, 13247, 12951, 2025, 9548, 2149, 2197, 9510, 6131, 5908, 803, 14342, 12995, 2064, 16024, 12992, 14332, 13202, 793, 2933, 11104, 379, 2056, 2055, 2922, 9506, 11344, 766, 4688, 6869, 6872, 3554, 2057, 2058, 8041, 9121, 9546, 7252, 7060, 9553, 5930, 5931, 5929, 2079, 5905, 9527, 4689, 3947, 1563, 9507, 4613, 794, 11775, 11271, 3951, 2062, 2027, 2060, 2061, 2070, 2032, 7154, 7153, 11779, 5932, 2976, 4613, 794, 1043, 4632, 4631, 4633, 4404, 2059, 12095, 2040, 2915, 4372, 2069, 11229, 7848, 15383, 10125, 2031, 10875, 12889, 12888, 321, 322, 9551, 9550, 3954, 2020, 2068, 4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552])

Main()