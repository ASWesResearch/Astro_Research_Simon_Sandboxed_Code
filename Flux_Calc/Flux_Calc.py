"""
This code is based on the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
"Exploring the High Energy Sky: A Survey of X-ray Sources in Galaxies Within 15 Megaparsecs of the Milky Way".
This thesis was submitted to the faculty of Wesleyan University in partial fulfillment of the requirements for the Degree of Bachelor of Arts
with Departmental Honors in Astronomy.

This version of the code is now modifed by Anthony Santini (Wesleyan 2019 Master's alumnus) for use in current research.


"""

import os
import sys
from ciao_contrib.runtool import *
import ciao_contrib.runtool as rt
import numpy as np
import glob
import re
import time
import logging
from multiprocessing import Pool
#Constants:
#Root_Path="/Volumes/"
Root_Path="/opt/"
#Reg_Path=Root_Path+"xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/Hybrid_Regions/"
Reg_Path="/Volumes/expansion/Hybrid_Regions/"
#Outroot_Path=root_Path+"xray/anthony/Simon_Sandboxed_Code/Flux_Calc/"
Outroot_Path="/Volumes/expansion/"
#Extracted_Path=Root_Path+"xray/anthony/Simon_Sandboxed_Code/Specextract/"
Extracted_Path="/Volumes/expansion/"
#Exception Classes
class FileNotFoundException(Exception):
    "Raised when a file is not found"
    def __init__(self, ObsID):
        self.ObsID = ObsID
        self.message = "ObsID "+str(self.ObsID)+" is Missing Evt2 File"
        super().__init__(self.message)

print(sys.version)
sys.path.append(Root_Path+"xray/anthony/Research_Git/")
from ObsID_From_CSV_Query import ObsID_From_CSV_Query
def Flux_Calc(obsid, Overlap_Only_Bool=True):
    print("ObsID: ", str(obsid))
    f = open("Srcflux_Times.txt", "a")
    f.write("ObsID,Source_Num,Time\n")
    #/Volumes/expansion/Source_Fluxes/10125/
    Error_Path="/Volumes/expansion/Source_Fluxes/"+str(obsid)+"/Flux_Calc_Error_Log.txt"
    directory = os.path.dirname(Error_Path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    Error_Log_File=open(Error_Path,"w")
    Error_List=[]
    with rt.new_pfiles_environment(ardlib=True):
        with new_tmpdir() as tmpdir:
            obsid = str(obsid)
            print("Beginning flux extraction for obsid " + obsid + "...")
            print("Getting name of evt2 file...")
            ##name = glob.glob("chandra_from_csc/" + obsid + "/primary/*_evt2.fits")
            #Evt2_Fpath_L = glob.glob(Root_Path+"xray/simon/all_chandra_observations/" + obsid + "/primary/*_evt2.fit*")
            Evt2_Fpath_L = glob.glob("/Volumes/expansion/ObsIDs/" + obsid + "/new/*_evt2.fit*")
            try:
                if(len(Evt2_Fpath_L)==0):
                    print("ObsID "+str(obsid)+" is Missing Evt2 File")
                    raise FileNotFoundException(obsid)
                Evt2_Fpath=Evt2_Fpath_L[0]
            except:
                print("ObsID "+str(obsid)+" is Missing Evt2 File")
                raise FileNotFoundException(obsid)
            print("Evt2_Fpath: ", Evt2_Fpath)
            ##print "Shortening name..."
            #Ant: This is done only for an exact directory on a computer not used anymore. It must be modified to create the shortened filename from an arbitary filepath
            ##name = name[0][-24:]
            name_L=Evt2_Fpath.split("/")
            print("name_L: ", name_L)
            name=name_L[len(name_L)-1]
            print(name)
            Hybrid_Filepath=Reg_Path+str(obsid)+"/"+str(obsid)+"_Nearest_Neighbor_Hybrid.reg"
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
            print("num: ", num)
            #print "Unlearning specextract..."
            #os.system("punlearn specextract")
            #print "Getting system path..."
            Hybrid_BG_Filepath=Reg_Path+str(obsid)+"/"+str(obsid)+"_Nearest_Neighbor_Hybrid_Background.reg"
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
            #specextract = rt.make_tool("specextract")
            dmcoords = rt.make_tool("dmcoords")
            srcflux = rt.make_tool("srcflux")
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
                Band_L.append(Band_L_Med)
                Bands_L.append(Band_L)
            print("Bands_L: ", Bands_L)
            #"0.5:7:2.5,ultrasoft"
            Bands_Str=""
            for Band_L in Bands_L:
                Cur_Band_Str=str(Band_L[0])+":"+str(Band_L[1])+":"+str(Band_L[2])+","
                Bands_Str=Bands_Str+Cur_Band_Str
            Bands_Str=Bands_Str+"csc"
            print("Bands_Str: ", Bands_Str)
            for i in range(0,num):
            #for i in range(57,62): #For Testing
                Hybrid_Region_Str=Hybrid_Region_Str_Reduced_L[i]
                print("Hybrid_Region_Str: ", Hybrid_Region_Str)
                Hybrid_Region_Str_L=re.split("[(),]", Hybrid_Region_Str)
                print("Hybrid_Region_Str_L: ", Hybrid_Region_Str_L)
                x_str=Hybrid_Region_Str_L[1]
                x=float(x_str)
                y_str=Hybrid_Region_Str_L[2]
                y=float(y_str)
                print("Got x and y values...")
                #print "punlearn dmcoords..."
                #os.system("punlearn dmcoords")
                #Infile_Str=Root_Path+"xray/simon/all_chandra_observations/" + str(obsid) + "/primary/" + str(name)
                Infile_Str="/Volumes/expansion/ObsIDs/" + str(obsid) + "/new/" + str(name)
                print("running dmcoords...")
                #print "dmcoords '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "' option=sky x=" + str(x) + " y=" + str(y) + " celfmt=deg"
                #os.system("dmcoords '/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "' option=sky x=" + str(x) + " y=" + str(y) + " celfmt=deg")
                dmcoords(infile=Infile_Str, x=x, y=y, option='sky', verbose=0, celfmt='deg')
                print("Getting Chip ID...")
                #chip = sp.check_output("pget dmcoords chip_id", shell=True)
                chip=dmcoords.chip_id
                chip = int(chip)
                print("Got Chip ID: " + str(chip))

                if chip == 5 or chip == 7:
                    inst = 's'
                else:
                    inst = 'i'
                #RA = sp.check_output("pget dmcoords ra", shell=True)
                RA=dmcoords.ra
                RA=float(RA)
                #Dec = sp.check_output("pget dmcoords dec", shell=True)
                Dec=dmcoords.dec
                Dec=float(Dec)
                print(inst)
                #Ant: Reducing the Hybrid_Region_Str to a region shape
                print("Extracting soft X-rays for source " + str(i+1) + " in obsid " + obsid + "...")
                Hybrid_Region_Str_Reduced_Split_L=re.split("[; ]",Hybrid_Region_Str)
                print("Hybrid_Region_Str_Reduced_Split_L: ", Hybrid_Region_Str_Reduced_Split_L)
                Hybrid_Region_Str_Reduced=Hybrid_Region_Str_Reduced_Split_L[1]
                print("Hybrid_Region_Str_Reduced: ", Hybrid_Region_Str_Reduced)
                #Ant: Parsing the parsing the Hybrid source background region file to in order to feed that input into srcflux
                j=i
                if(j%2!=0):
                    j=j+1
                Hybrid_BG_Region_Outer_R_Str=Hybrid_BG_Region_Str_Reduced_L[j]
                Hybrid_BG_Region_Outer_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Outer_R_Str)
                Hybrid_BG_Region_Outer_R_Str_Reduced=Hybrid_BG_Region_Outer_R_Str_Reduced_L[1]
                print("Hybrid_BG_Region_Outer_R_Str_Reduced: ", Hybrid_BG_Region_Outer_R_Str_Reduced)
                Hybrid_BG_Region_Inner_R_Str=Hybrid_BG_Region_Str_Reduced_L[j+1]
                Hybrid_BG_Region_Inner_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Inner_R_Str)
                Hybrid_BG_Region_Inner_R_Str_Reduced=Hybrid_BG_Region_Inner_R_Str_Reduced_L[1]
                print("Hybrid_BG_Region_Inner_R_Str_Reduced: ", Hybrid_BG_Region_Inner_R_Str_Reduced)
                Hybrid_BG_Region=Hybrid_BG_Region_Outer_R_Str_Reduced+Hybrid_BG_Region_Inner_R_Str_Reduced
                #For overlapping sources
                Hybrid_Overlap_Corrected_BG_Filepath=Reg_Path+str(obsid)+"/Individual_Source_Regions/Source_"+str(i+1)+"_ObsID_"+str(obsid)+"_Nearest_Neighbor_Hybrid_Overlap_Corrected_Background.reg"
                Hybrid_Overlap_Corrected_BG_File=open(Hybrid_Overlap_Corrected_BG_Filepath)
                Hybrid_Overlap_Corrected_BG_Str=Hybrid_Overlap_Corrected_BG_File.read()
                Hybrid_Overlap_Corrected_BG_Str_L=Hybrid_Overlap_Corrected_BG_Str.split("\n")
                print("Hybrid_Overlap_Corrected_BG_Str_L: ", Hybrid_Overlap_Corrected_BG_Str_L)
                print("len(Hybrid_Overlap_Corrected_BG_Str_L): ", len(Hybrid_Overlap_Corrected_BG_Str_L))
                if(len(Hybrid_Overlap_Corrected_BG_Str_L)>5):
                    Overlap_Bool=True
                else:
                    Overlap_Bool=False
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
                print("Running srcflux...")
                print("Source " + str(i+1) + " of " + str(num) + " for obsid " + str(obsid) + "...")
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
                """
                os.system("specextract \
                infile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_Region_Str_Reduced+"]\" \
                bkgfile=\"/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_BG_Region+"]\" \
                outroot=\"/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/" + str(obsid) +"/"+ "extracted_spectra" + obsid + "_" + str(i+1) + "' \
                correctpsf=yes \
                weight=no \
                clobber=yes \
                verbose=1")
                """

                """
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
                """
                """
                srcflux  infile pos outroot [bands] [srcreg] [bkgreg] [bkgresp]
                [psfmethod] [psffile] [conf] [binsize] [rmffile] [arffile] [model]
                [paramvals] [absmodel] [absparams] [abund] [fovfile] [asolfile]
                [mskfile] [bpixfile] [dtffile] [ecffile] [parallel] [nproc] [tmpdir]
                [clobber] [verbose]
                """
                """
                #Infile_Str="/Volumes/xray/simon/all_chandra_observations/" + str(obsid) + "/primary/" + str(name)
                ##Bkgfile_Str="/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_BG_Region+"]"
                Outroot_Str=Root_Path+"xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes/" + str(obsid) +"/"+ "Source_Flux_" + str(obsid) + "_" + str(i+1)
                Outroot_Str_Test=Root_Path+"xray/anthony/Simon_Sandboxed_Code/Flux_Calc/Source_Fluxes_Test/" + str(obsid) +"/"+ "Source_Flux_" + str(obsid) + "_" + str(i+1)
                #specextract(infile=Infile_Str, bkgfile=Bkgfile_Str, outroot=Outroot_Str, correctpsf="yes", weight="no", clobber="yes", tmpdir=tmpdir, verbose=1)
                Rmffile_Str=Root_Path+"xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".rmf"
                Arffile_Str=Root_Path+"xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".arf"
                """
                #Infile_Str="/Volumes/xray/simon/all_chandra_observations/" + str(obsid) + "/primary/" + str(name)
                ##Bkgfile_Str="/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_BG_Region+"]"
                Outroot_Str=Outroot_Path+"Source_Fluxes/" + str(obsid) +"/"+ "Source_Flux_" + str(obsid) + "_" + str(i+1)
                Outroot_Str_Test=Outroot_Path+"/Source_Fluxes_Test/" + str(obsid) +"/"+ "Source_Flux_" + str(obsid) + "_" + str(i+1)
                #specextract(infile=Infile_Str, bkgfile=Bkgfile_Str, outroot=Outroot_Str, correctpsf="yes", weight="no", clobber="yes", tmpdir=tmpdir, verbose=1)
                Rmffile_Str=Extracted_Path+"extracted_spectra/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".rmf"
                Arffile_Str=Extracted_Path+"extracted_spectra/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".arf"
                #/Volumes/expansion/extracted_spectra_Test/
                Rmffile_Str_Test=Extracted_Path+"extracted_spectra_Test/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".rmf"
                Arffile_Str_Test=Extracted_Path+"extracted_spectra_Test/"+ obsid +"/extracted_spectra_"+ obsid +"_"+str(i+1)+".arf"
                #Bkgfile_Str_Overlap="/Volumes/expansion/ObsIDs/" + obsid + "/new/" + str(name) + "[sky=region("+Hybrid_Overlap_Corrected_BG_Filepath+")]"
                Bkgfile_Str_Overlap="region("+Hybrid_Overlap_Corrected_BG_Filepath+")"
                print("Bkgfile_Str_Overlap: ", Bkgfile_Str_Overlap)

                print(type(Dec))
                print("Dec<0: ", Dec<0)
                if(Dec>0):
                    Pos_Str=str(RA)+" +"+str(Dec)
                elif(Dec<0):
                    Pos_Str=str(RA)+" "+str(Dec)
                #Pos_Str=str(RA)+","+str(Dec)
                print("Pos_Str: ", Pos_Str)
                #Ant: Note I need to fix the bands parameter. The string imput is incorrect and using the incorrect syntax ! ! !, Ant Note 2: This has now been corrected !
                tic = time.time()
                print("Overlap_Bool: ", Overlap_Bool)
                try:
                    if(Overlap_Only_Bool and Overlap_Bool):
                        print("Overlap Detected!")
                        #srcflux(infile=Infile_Str, pos=Pos_Str, outroot=Outroot_Str_Test, bands=Bands_Str, srcreg=Hybrid_Region_Str_Reduced, bkgreg=Bkgfile_Str_Overlap, psfmethod="ideal", rmffile=Rmffile_Str_Test, arffile=Arffile_Str_Test, clobber="yes", tmpdir=tmpdir, verbose=5)
                        #srcflux(infile=Infile_Str, pos=Pos_Str, outroot=Outroot_Str, bands=Bands_Str, srcreg=Hybrid_Region_Str_Reduced, bkgreg=Bkgfile_Str_Overlap, psfmethod="ideal", rmffile=Rmffile_Str, arffile=Arffile_Str, clobber="yes", tmpdir=tmpdir, verbose=5)
                        srcflux(infile=Infile_Str, pos=Pos_Str, outroot=Outroot_Str, bands=Bands_Str, srcreg=Hybrid_Region_Str_Reduced, bkgreg=Bkgfile_Str_Overlap, psfmethod="ideal", rmffile=Rmffile_Str, arffile=Arffile_Str, clobber="no", tmpdir=tmpdir, verbose=1)
                    else:
                        print("No Clobber")
                        #srcflux(infile=Infile_Str, pos=Pos_Str, outroot=Outroot_Str_Test, bands=Bands_Str, srcreg=Hybrid_Region_Str_Reduced, bkgreg=Hybrid_BG_Region, psfmethod="ideal", rmffile=Rmffile_Str_Test, arffile=Arffile_Str_Test, clobber="no", tmpdir=tmpdir, verbose=5)
                        srcflux(infile=Infile_Str, pos=Pos_Str, outroot=Outroot_Str, bands=Bands_Str, srcreg=Hybrid_Region_Str_Reduced, bkgreg=Hybrid_BG_Region, psfmethod="ideal", rmffile=Rmffile_Str, arffile=Arffile_Str, clobber="no", tmpdir=tmpdir, verbose=1)
                except Exception as Argument:
                    Error_List.append([obsid,i+1])
                    Error_Log_File.write(str(obsid)+" "+"Srcflux Error:\n"+str(Argument)+"\n")
                toc = time.time()
                #print(f"Time for srcflux: {toc - tic:0.4f} seconds")
                Time=toc - tic
                print("Time for srcflux for ObsID "+str(obsid)+" Source " + str(i+1)+":", Time)
                f.write(str(obsid)+","+str(i+1)+","+str(Time)+"\n")
    print("Error_List: ", Error_List)
    Error_Log_File.write("\n\n\nFail_List: "+str(Error_List))
    Error_Log_File.close()
    f.close()
    Hybrid_File.close()
    Hybrid_BG_File.close()
    Hybrid_Overlap_Corrected_BG_File.close()



##if __name__ == '__main__':
def Driver(obsid_L):
    """
    for i in range(0,len(obsid_L)):
        try:
            #print "Attempting specextract..."
            spec(obsid_L[i])
        except OSError:
            print "No traced sources in obsid " + str(obsid[i])
            print "specextract failed for obsid " + str(obsid[i])
    """
    """
    for i in range(0,len(obsid_L)):
        Flux_Calc(obsid_L[i])
    """
    #"""
    if __name__ == '__main__':
        P = Pool()
        P.map(Flux_Calc, obsid_L)
    #"""
    #when_done("specextract") #Ant: Disabled to not accidently email Simon

def Main():
    #Driver([10125])
    #Driver([11775]) #NGC_1300 to test negitive declination
    #Driver([2031])
    ##Driver([6096, 1971, 1972, 768, 952, 11674, 13255, 13253, 13246, 12952, 12953, 13247, 12951, 2025, 9548, 2149, 2197, 9510, 6131, 5908, 803, 14342, 12995, 2064, 16024, 12992, 14332, 13202, 793, 2933, 11104, 379, 2056, 2055, 2922, 9506, 11344, 766, 4688, 6869, 6872, 3554, 2057, 2058, 8041, 9121, 9546, 7252, 7060, 9553, 5930, 5931, 5929, 2079, 5905, 9527, 4689, 3947, 1563, 9507, 4613, 794, 11775, 11271, 3951, 2062, 2027, 2060, 2061, 2070, 2032, 7154, 7153, 11779, 5932, 2976, 4613, 794, 1043, 4632, 4631, 4633, 4404, 2059, 12095, 2040, 2915, 4372, 2069, 11229, 7848, 15383, 10125, 2031, 10875, 12889, 12888, 321, 322, 9551, 9550, 3954, 2020, 2068, 4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552])
    ##ObsID_L=ObsID_From_CSV_Query.Read_ObsIDs(Remove_Unarchived=True)
    ##Driver(ObsID_L)
    #Driver([14676, 349, 350, 353, 354, 361, 16745, 378, 379, 380, 12668, 383, 384, 388, 389, 390, 391, 392, 393, 394, 395, 400, 402, 404, 405, 407, 12696, 409, 410, 411, 24981, 413, 414, 24986, 18875, 4555, 4556, 4557, 4558, 12748, 14795, 14801, 10722, 10723, 10724, 10725, 10726, 20965, 20966, 20992, 20993, 4613, 20997, 20998, 20999, 21000, 21001, 21003, 4627, 4628, 4629, 4630, 23075, 23076, 21036, 14896, 14902, 14912, 6727, 16969, 4688, 4689, 4690, 16978, 4692, 4693, 4694, 16983, 4696, 4697, 21077, 21082, 16991, 16994, 16995, 16996, 16997, 23140, 23141, 17000, 18440, 17003, 17007, 10868, 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 2686, 4734, 4735, 4736, 4737, 6781, 6782, 4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 10925, 23216, 23217, 12978, 23218, 23219, 12981, 23220, 23223, 12992, 12993, 12994, 12995, 12996, 13018, 735, 23266, 21230, 17155, 782, 784, 790, 792, 793, 794, 795, 11032, 797, 11033, 11034, 17180, 808, 15149, 2879, 2885, 11080, 11081, 11082, 11083, 11084, 11085, 11086, 15190, 864, 11104, 15200, 19297, 2916, 2917, 870, 871, 872, 2918, 2919, 19304, 21350, 2925, 21351, 21352, 882, 2933, 2934, 2949, 2950, 21384, 19339, 19344, 19345, 13202, 19346, 7060, 19348, 19350, 19351, 19354, 7069, 19357, 2976, 7073, 9122, 2978, 7074, 7075, 7076, 934, 9120, 9121, 7082, 7083, 7084, 942, 7086, 7087, 19374, 7090, 7091, 23472, 7093, 23473, 7095, 7096, 13241, 7098, 19386, 19387, 7101, 15294, 7103, 7104, 7105, 962, 963, 3012, 7106, 13248, 7111, 15295, 969, 7113, 7115, 7116, 19397, 7118, 19403, 7120, 7121, 19407, 7123, 7124, 19411, 19414, 7127, 19416, 19417, 7132, 19421, 7134, 19422, 21471, 21472, 21473, 21474, 19428, 15333, 21479, 7146, 7147, 19437, 7150, 7152, 7153, 7154, 13303, 13304, 11260, 11268, 23559, 11272, 11273, 23561, 23564, 15382, 15384, 11289, 15386, 15387, 11295, 19497, 21545, 11309, 11311, 11317, 17461, 17462, 9278, 17471, 17472, 19521, 19522, 19524, 5197, 7252, 25689, 13439, 21639, 15496, 21640, 17547, 19363, 21647, 21648, 21649, 17569, 17570, 5283, 17571, 17578, 5296, 5297, 5300, 5301, 5302, 5309, 15553, 21698, 21699, 7369, 5322, 5323, 15572, 15574, 5337, 5338, 5339, 5340, 15579, 15582, 15587, 15588, 15589, 15594, 15603, 3325, 15616, 17678, 1302, 15646, 19392, 19747, 19393, 19748, 19394, 9532, 9533, 9534, 9535, 9536, 9537, 9538, 9539, 9540, 9541, 9542, 9543, 9545, 9546, 9547, 9548, 9549, 9550, 9551, 9552, 9553, 23474, 21853, 23475, 9570, 23476, 23477, 23478, 23479, 13686, 23480, 23481, 23482, 23483, 23484, 15756, 23485, 15760, 23486, 23487, 15771, 23488, 13728, 23489, 23490, 23491, 23492, 23493, 10875, 15803, 23494, 23495, 13765, 23496, 23497, 3550, 3551, 13791, 17890, 17891, 13796, 11761, 5619, 13812, 13813, 13814, 13815, 13816, 13817, 13819, 13820, 13821, 13822, 13829, 11782, 13830, 13831, 13832, 11786, 5644, 19981, 19982, 11800, 1564, 1578, 1579, 1586, 1587, 11846, 11847, 9805, 1618, 1621, 1622, 1624, 1633, 1634, 1635, 1636, 1637, 1638, 1640, 7797, 7798, 7799, 7800, 14984, 18047, 16000, 16001, 14985, 16002, 16003, 16005, 18048, 18053, 18054, 18062, 18063, 18064, 18065, 18066, 18067, 18068, 9877, 18069, 16023, 16024, 18070, 18071, 9883, 16028, 16029, 18072, 18073, 16032, 16033, 7850, 22189, 7858, 22194, 7863, 17032, 14017, 14018, 16068, 16069, 3786, 3787, 3788, 7885, 16121, 16122, 5905, 5911, 5929, 5930, 5931, 10025, 10026, 10027, 5935, 5936, 5937, 5938, 5939, 5940, 5941, 5942, 5943, 5944, 5945, 5946, 5947, 5948, 5949, 12095, 3925, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 3937, 3938, 3939, 3940, 3941, 3942, 3943, 22372, 22375, 16234, 3949, 3950, 20333, 3953, 3954, 8050, 8052, 8053, 20343, 8058, 12155, 12156, 3965, 20353, 16260, 16261, 16262, 20356, 10125, 16276, 16277, 8086, 14230, 14231, 8091, 8098, 18340, 18341, 18342, 18343, 4010, 4016, 4017, 18352, 4019, 8125, 8126, 12238, 12239, 6096, 6097, 22478, 22479, 22480, 22481, 22482, 2014, 6114, 6115, 6118, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2039, 2040, 14332, 8190]) #Note: This may not be the full sample.
    #Driver([4613])
    #Driver([969])
    ##print("ObsID_L:\n", ObsID_L)
    ##ObsID_Mando_L=[10274, 2014, 2065, 2917, 353, 3949, 4628, 4735, 5283, 7087, 11080, 2025, 2075, 2918, 354, 3950, 4629, 4736, 5296, 7252, 12473, 2026, 2076, 2919, 3550, 3953, 4693, 4737, 5297, 735, 13819, 2027, 21474, 2925, 3551, 3954, 4694, 4743, 5300, 782, 14801, 2028, 2148, 2933, 361, 3965, 4696, 4744, 5301, 784, 15760, 2029, 2197, 2934, 3786, 4010, 4697, 4745, 5302, 790, 1633, 2030, 2198, 2949, 3787, 4016, 4725, 4746, 5309, 792, 1634, 2031, 2255, 2950, 3788, 4017, 4726, 4747, 5322, 793, 1635, 2032, 2260, 2976, 3925, 4019, 4727, 4748, 5323, 9278, 1636, 20333, 23220, 2978, 3937, 4176, 4728, 4749, 5945, 1637, 2039, 2340, 3012, 3938, 4555, 4729, 4750, 5946, 1638, 2040, 23499, 316, 3939, 4556, 4730, 4751, 5947, 1640, 2057, 2686, 318, 3940, 4557, 4731, 4752, 5948, 16978, 2058, 2879, 3325, 3941, 4558, 4732, 4753, 5949, 18063, 2059, 2885, 349, 3942, 4613, 4733, 4754, 6096, 19350, 2064, 2916, 350, 3943, 4627, 4734, 5197, 6097, 378, 379, 380, 383, 384, 388, 389, 390, 391, 392, 393, 394, 395, 400, 402, 404, 405, 407, 409, 410, 411, 413, 414, 794, 795, 797, 808, 864, 870, 871, 872, 882, 934, 942, 962, 963, 969, 1302, 1564, 1578, 1579, 1586, 1587, 1618, 1621, 1622, 1624, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 4630, 4688, 4689, 4690, 4692, 5337, 5338, 5339, 5340, 5619, 5644, 5905, 5911, 5929, 5930, 5931, 5935, 5936, 5937, 5938, 5939, 5940, 5941, 5942, 5943, 5944, 6114, 6115, 6118, 6152, 6169, 6170, 6175, 6184, 6185, 6361, 6727, 6781, 6782, 7060, 7069, 7073, 7074, 7075, 7076, 7082, 7083, 7084, 7086, 7090, 7091, 7093, 7095, 7096, 7098, 7101, 7103, 7104, 7105, 7106, 7111, 7113, 7115, 7116, 7118, 7120, 7121, 7123, 7124, 7127, 7132, 7134, 7146, 7147, 7150, 7152, 7153, 7154, 7369, 7797, 7798, 7799, 7800, 7850, 7858, 7863, 7885, 8050, 8052, 8053, 8058, 8086, 8091, 8098, 8125, 8126, 8190, 8197, 8198, 8458, 8464, 8465, 8489, 8490, 9120, 9121, 9122, 9532, 9533, 9534, 9535, 9536, 9537, 9538, 9539, 9540, 9541, 9542, 9543, 9545, 9546, 9547, 9548, 9549, 9550, 9551, 9552, 9553, 9570]
    #Remainder of Mando_L
    ##ObsID_Mando_L_Remainder=[383,384,388,389,390,391,392,393,394,395,400,402,963,1302,1564,1578,1579,1586,1587,4692,5337,5338,5339,5340,5619,5644,5905,5911,5929,5930,5931,5935,7799,7800,7850,7858,7863,7885]
    #ObsID_Mando_L_Remainder=[383,384,388,389,390,391,392,393,394,395,402,1302,1564,1579,1586,1587,4692,5337,5338,5339,5340,5619,5644,5905,5911,5929,5930,5931,5935,7799,7800,7850,7858,7863,7885] #Removed ObsIDs 400,963,1578 due to missing files.
    #Vetinari_L_Sub_1=[9805, 9877, 9883, 10025, 10026, 10027, 10125, 10275, 10276, 10277, 10278, 10279, 10280, 10281, 10289, 10290, 10291, 10292, 10293, 10542, 10543, 10544, 10545, 10559, 10560, 10722, 10723, 10724, 10725, 10726]
    #Vetinari_L_Sub_2=[10868, 10875, 10925, 11032, 11033, 11034, 11081, 11082, 11083, 11084, 11085, 11086, 11104, 11260, 11268, 11272, 11273, 11289, 11295, 11309, 11311, 11317, 11761, 11782, 11786, 11800, 11846, 11847, 12095, 12155]
    #Vetinari_L_Sub_3=[12156, 12238, 12239, 12301, 12437, 12562, 12668, 12696, 12748, 12978, 12981, 12992, 12993, 12994, 12995, 12996, 13018, 13202, 13241, 13248, 13303, 13304, 13439, 13686, 13728, 13765, 13791, 13796, 13812, 13813]
    #Vetinari_L_Sub_4=[13814, 13815, 13816, 13817, 13820, 13821, 13822, 13829, 13830, 13831, 13832, 14017, 14018, 14230, 14231, 14332, 14341, 14342, 14349, 14350, 14351, 14376, 14378, 14383, 14384, 14412, 14419, 14437, 14442, 14471]
    #Vetinari_L_Sub_5=[14675, 14676, 14795, 14896, 14902, 14912, 14984, 14985, 15149, 15190, 15200, 15294, 15295, 15333, 15382, 15384, 15386, 15387, 15496, 15553, 15572, 15574, 15579, 15582, 15587, 15588, 15589, 15594, 15603, 15616]
    #Vetinari_L_Sub_6=[15646, 15756, 15771, 15803, 16000, 16001, 16002, 16003, 16005, 16023, 16024, 16028, 16029, 16032, 16033, 16068, 16069, 16121, 16122, 16234, 16260, 16261, 16262, 16276, 16277, 16484, 16485, 16556, 16580, 16745]
    #Vetinari_L_Sub_7=[16969, 16983, 16991, 16994, 16995, 16996, 16997, 17000, 17003, 17007, 17032, 17155, 17180, 17461, 17462, 17471, 17472, 17547, 17569, 17570, 17571, 17578, 17678, 17890, 17891, 18047, 18048, 18053, 18054, 18062]
    #Vetinari_L_Sub_8=[18064, 18065, 18066, 18067, 18068, 18069, 18070, 18071, 18072, 18073, 18340, 18341, 18342, 18343, 18352, 18440, 18454, 18455, 18461, 18462, 18760, 18875, 19297, 19304, 19339, 19344, 19345, 19346, 19348, 19351]
    #Vetinari_L_Sub_9=[19354, 19357, 19363, 19374, 19386, 19387, 19392, 19393, 19394, 19397, 19403, 19407, 19411, 19414, 19416, 19417, 19421, 19422, 19428, 19437, 19497, 19521, 19522, 19524, 19747, 19748, 19981, 19982, 20343, 20353]
    #Vetinari_L_Sub_11=[21640, 21647, 21648, 21649, 21698, 21699, 21853, 22189, 22194, 22372, 22375, 22478, 22479, 22480, 22481, 22482, 22714, 22715, 23075, 23076, 23140, 23141, 23216, 23217, 23218, 23219, 23223, 23266, 23472, 23473]
    #Vetinari_L_Sub_12=[23474, 23475, 23476, 23477, 23478, 23479, 23480, 23481, 23482, 23483, 23484, 23485, 23486, 23487, 23488, 23489, 23490, 23491, 23492, 23493, 23494, 23495, 23496, 23497, 23498, 23500, 23501, 23559, 23561, 23564]
    #Vetinari_L_Sub_12=[23474, 23475, 23476, 23477, 23478, 23479, 23480, 23481, 23482, 23483, 23484, 23485, 23486, 23487, 23488, 23489, 23490, 23491, 23492, 23493, 23494, 23495, 23496, 23497, 23498, 23500, 23501, 23559, 23561, 23564]
    Vetinari_L_Sub_12=[23474, 23475, 23476, 23477, 23478, 23479, 23480, 23481, 23482, 23483, 23484, 23485, 23486, 23487, 23488, 23489, 23490, 23491, 23492, 23493, 23494, 23495, 23496, 23497, 23498, 23559, 23561, 23564] #Unarchived ObsIDs removed.
    ##Driver(ObsID_Mando_L)
    #Driver(ObsID_Mando_L_Remainder)
    #Driver(Vetinari_L_Sub_1)
    #Driver(Vetinari_L_Sub_2)
    #Driver(Vetinari_L_Sub_4)
    #Driver(Vetinari_L_Sub_5)
    #Driver(Vetinari_L_Sub_6)
    #Driver(Vetinari_L_Sub_7)
    #Driver(Vetinari_L_Sub_8)
    #Driver(Vetinari_L_Sub_9)
    #Driver(Vetinari_L_Sub_11)
    Driver(Vetinari_L_Sub_12)
Main()
