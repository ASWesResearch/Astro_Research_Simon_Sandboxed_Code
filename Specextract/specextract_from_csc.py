"""
This code is the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
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
from multiprocessing import Pool
#from done_email import when_done
path="/opt/xray/anthony/Research_Git/"
sys.path.append(os.path.abspath(path))
from ObsID_From_CSV_Query import ObsID_From_CSV_Query #Ant: This should be changed as to not depend on the Research_Git directory
#Constants:
#Root_Path="/Volumes/"
Root_Path="/opt/"
#Reg_Path="/Volumes/xray/anthony/Research_Git/Nearest_Raytraced_Neighbor_Calc/"
Reg_Path="/Volumes/expansion/"
#Outpath="/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/"
Outpath="/Volumes/expansion/"
def spec(obsid, Clobber_Bool=False):
    Error_List=[]
    with rt.new_pfiles_environment(ardlib=True):
        with new_tmpdir() as tmpdir:
            """
            # tmpdir has been created and can be used in this block
            bmap = make_tool("create_bkg_map")
            bmap.tmpdir = tmpdir
            ...
            bmap(...)
            """
            """
            print "Beginning specextract for obsid " + str(obsid)
            obsid = str(obsid)
            print "Getting name of evt2 file..."
            name = glob.glob("chandra_from_csc/" + obsid + "/primary/*_evt2.fits")
            print name
            print "Shortening name..."
            name = name[0][-24:]
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
            """
            obsid = str(obsid)
            #/Volumes/expansion/extracted_spectra/10125/
            if(Clobber_Bool==False):
                Output_Path="/Volumes/expansion/extracted_spectra/"+str(obsid)+"/"
                Output_Path_L=glob.glob(Output_Path)
                if(len(Output_Path_L)>0):
                    print(str(obsid)+" already exists and Clobber is set to False!")
                    return
            print("Beginning flux extraction for obsid " + obsid + "...")
            print("Getting name of evt2 file...")
            ##name = glob.glob("chandra_from_csc/" + obsid + "/primary/*_evt2.fits")
            #Evt2_Fpath_L = glob.glob("/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/*_evt2.fit*")
            Evt2_Fpath_L = glob.glob("/Volumes/expansion/ObsIDs/" + obsid + "/new/*evt2.fit*")
            if(len(Evt2_Fpath_L)==0):
                print(str(obsid)+" has no reprocessed Evt2 file!")
                Evt2_Fpath_L = glob.glob("/Volumes/expansion/ObsIDs/" + obsid + "/*evt2.fit*")
            if(len(Evt2_Fpath_L)==0):
                print(str(obsid)+" HAS NO MATCHING ObsIDs!!!!")
                return
            Evt2_Fpath=Evt2_Fpath_L[0]
            print("Evt2_Fpath: ", Evt2_Fpath)
            ##print "Shortening name..."
            #Ant: This is done only for an exact directory on a computer not used anymore. It must be modified to create the shortened filename from an arbitary filepath
            ##name = name[0][-24:]
            name_L=Evt2_Fpath.split("/")
            print("name_L: ", name_L)
            name=name_L[len(name_L)-1]
            print(name)
            Hybrid_Filepath=Reg_Path+"Hybrid_Regions/"+str(obsid)+"/"+str(obsid)+"_Nearest_Neighbor_Hybrid.reg"
            try:
                Hybrid_File=open(Hybrid_Filepath)
            except:
                print(str(obsid)+" has no source region file!")
                return
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
            ##print "Unlearning specextract..."
            ##os.system("punlearn specextract")
            #print "Getting system path..."
            Hybrid_BG_Filepath=Reg_Path+"Hybrid_Regions/"+str(obsid)+"/"+str(obsid)+"_Nearest_Neighbor_Hybrid_Background.reg"
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
            specextract = rt.make_tool("specextract")
            print("specextract Tool Object Before:\n", specextract)
            #specextract.tmpdir = tmpdir
            print("tempdir: " , tmpdir)
            specextract.tmpdir = tmpdir
            print("specextract Tool Object After:\n", specextract)
            #for i in range(1,num+1):
            for i in range(0,num):
                Hybrid_Region_Str=Hybrid_Region_Str_Reduced_L[i]
                print("Hybrid_Region_Str: ", Hybrid_Region_Str)
                Hybrid_Region_Str_L=re.split("[(),]", Hybrid_Region_Str)
                print("Hybrid_Region_Str_L: ", Hybrid_Region_Str_L)
                """
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
                """
                #Ant: Reducing the Hybrid_Region_Str to a region shape
                print("Extracting soft X-rays for source " + str(i+1) + " in obsid " + obsid + "...")
                Hybrid_Region_Str_Reduced_Split_L=re.split("[; ]",Hybrid_Region_Str)
                print("Hybrid_Region_Str_Reduced_Split_L: ", Hybrid_Region_Str_Reduced_Split_L)
                Hybrid_Region_Str_Reduced=Hybrid_Region_Str_Reduced_Split_L[1]
                print("Hybrid_Region_Str_Reduced: ", Hybrid_Region_Str_Reduced)
                #Ant: Parsing the parsing the Hybrid source background region file to in order to feed that input into specextract
                print("i: ", i)
                j=i
                print("j Before: ", j)
                """
                if(j%2!=0):
                    #j=j+1 #Ant: I lost a week to this bug...
                    j=2*j
                """
                j=2*j
                print("j After: ", j)
                Hybrid_BG_Region_Outer_R_Str=Hybrid_BG_Region_Str_Reduced_L[j]
                Hybrid_BG_Region_Outer_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Outer_R_Str)
                Hybrid_BG_Region_Outer_R_Str_Reduced=Hybrid_BG_Region_Outer_R_Str_Reduced_L[1]
                print("Hybrid_BG_Region_Outer_R_Str_Reduced: ", Hybrid_BG_Region_Outer_R_Str_Reduced)
                Hybrid_BG_Region_Inner_R_Str=Hybrid_BG_Region_Str_Reduced_L[j+1]
                Hybrid_BG_Region_Inner_R_Str_Reduced_L=re.split("[; ]",Hybrid_BG_Region_Inner_R_Str)
                Hybrid_BG_Region_Inner_R_Str_Reduced=Hybrid_BG_Region_Inner_R_Str_Reduced_L[1]
                print("Hybrid_BG_Region_Inner_R_Str_Reduced: ", Hybrid_BG_Region_Inner_R_Str_Reduced)
                Hybrid_BG_Region=Hybrid_BG_Region_Outer_R_Str_Reduced+Hybrid_BG_Region_Inner_R_Str_Reduced
                print("Hybrid_BG_Region: ", Hybrid_BG_Region)
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
                print("Running specextract...")
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
                os.system("specextract \
                infile='/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_Region_Str_Reduced+"]' \
                bkgfile='/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_BG_Region+"]' \
                outroot='/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/" + str(obsid) +"/"+ "extracted_spectra_" + obsid + "_" + str(i+1) + "' \
                correctpsf=yes \
                weight=no \
                clobber=yes \
                verbose=1")
                """

                #Ant Note: dmcoords(infile=str(evtfpath),chipx=BG_X, chipy=BG_Y, chip_id=Chip_ID_Test, option='chip', verbose=0)
                """
                Infile_Str="/Volumes/xray/simon/all_chandra_observations/" + str(obsid) + "/primary/" + str(name) + "[sky="+Hybrid_Region_Str_Reduced+"]"
                Bkgfile_Str="/Volumes/xray/simon/all_chandra_observations/" + obsid + "/primary/" + str(name) + "[sky="+Hybrid_BG_Region+"]"
                Outroot_Str="/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra/" + str(obsid) +"/"+ "extracted_spectra_" + str(obsid) + "_" + str(i+1)
                Outroot_Str_Test="/Volumes/xray/anthony/Simon_Sandboxed_Code/Specextract/extracted_spectra_Test/" + str(obsid) +"/"+ "extracted_spectra_" + str(obsid) + "_" + str(i+1)
                """
                #specextract(infile=Infile_Str, bkgfile=Bkgfile_Str, outroot=Outroot_Str, correctpsf="yes", weight="no", clobber="yes", verbose=1)
                #Ant Note: ah = rt.make_tool("asphist")
                Infile_Str="/Volumes/expansion/ObsIDs/" + str(obsid) + "/new/" + str(name) + "[sky="+Hybrid_Region_Str_Reduced+"]"
                Bkgfile_Str="/Volumes/expansion/ObsIDs/" + obsid + "/new/" + str(name) + "[sky="+Hybrid_BG_Region+"]"
                Outroot_Str=Outpath+"extracted_spectra/" + str(obsid) +"/"+ "extracted_spectra_" + str(obsid) + "_" + str(i+1)
                Outroot_Str_Test=Outpath+"extracted_spectra_Test/" + str(obsid) +"/"+ "extracted_spectra_" + str(obsid) + "_" + str(i+1)
                """
                specextract = rt.make_tool("specextract")
                print "specextract Tool Object: ", specextract
                #specextract.tmpdir = tmpdir
                print "tempdir: " , tmpdir
                """
                print("Infile_Str: ", Infile_Str)
                print("Bkgfile_St: ", Bkgfile_Str)
                #print "Outroot_Str_Test: ", Outroot_Str_Test
                print("Outroot_Str: ", Outroot_Str)
                print("tempdir In: ", tmpdir)
                #Temp_Files = glob.glob(tmpdir+"| grep 'asphist'")
                #Temp_Files = glob.glob(tmpdir+"/*")
                #print "Temp_Files:\n", Temp_Files
                """
                with open(tmpdir, 'r') as f:
                    read_data = f.read()
                    f.closed
                """
                try:
                    #specextract(infile=Infile_Str, bkgfile=Bkgfile_Str, outroot=Outroot_Str_Test, correctpsf="yes", weight="no", clobber="yes", tmpdir=tmpdir, verbose=1)
                    if(Clobber_Bool):
                        specextract(infile=Infile_Str, bkgfile=Bkgfile_Str, outroot=Outroot_Str, correctpsf="yes", weight="no", clobber="yes", tmpdir=tmpdir, verbose=1)
                    else:
                        specextract(infile=Infile_Str, bkgfile=Bkgfile_Str, outroot=Outroot_Str, correctpsf="yes", weight="no", clobber="no", tmpdir=tmpdir, verbose=1)
                except:
                    Cur_Error_L=[obsid, Infile_Str, Bkgfile_Str, Outroot_Str_Test, tmpdir]
                    Error_List.append(Cur_Error_L)

    print("Error_List: ", Error_List)


##if __name__ == '__main__':
def Driver(obsid_L):
    """
    print "Reading in obsids... "
    obsid_path = glob.glob("chandra_from_csc/*")

    obsid = []
    for i in range(0,len(obsid_path)):
        #print obsid_path[i]
        obsid.append(obsid_path[i][17:])
    print "...done!"

    #print obsid
    #print len(obsid)
    obsid = clean_results(obsid)

    #print len(obsid)
    #print obsid

    print "All obsid names processed."
    """

    """
    processes=[]
    for _ in range(0,10):
        p=multiprocessing.Process(target=do_somthing, args=[1.5])
        p.start()
        processes.appned(p)
    for process in processes:
        process.join()
    """
    """
    from multiprocessing import Pool

    def f(x):
        return x*x

    if __name__ == '__main__':
        P = Pool(5)
        print(P.map(f, [1, 2, 3]))
    """

    """
    for i in range(0,len(obsid_L)):
        try:
            #print "Attempting specextract..."
            spec(obsid_L[i])
        except OSError:
            print "No traced sources in obsid " + str(obsid[i])
            print "specextract failed for obsid " + str(obsid[i])
    """
    #Ant: Multiprocessing to make it stop
    if __name__ == '__main__':
        P = Pool()
        P.map(spec, obsid_L)

    #when_done("specextract") #Ant: Disabled to not accidently email Simon

def Main():
    #Driver([10125])
    #Driver([6096, 1971, 1972, 768, 952, 11674, 13255, 13253, 13246, 12952, 12953, 13247, 12951, 2025, 9548, 2149, 2197, 9510, 6131, 5908, 803, 14342, 12995, 2064, 16024, 12992, 14332, 13202, 793, 2933, 11104, 379, 2056, 2055, 2922, 9506, 11344, 766, 4688, 6869, 6872, 3554, 2057, 2058, 8041, 9121, 9546, 7252, 7060, 9553, 5930, 5931, 5929, 2079, 5905, 9527, 4689, 3947, 1563, 9507, 4613, 794, 11775, 11271, 3951, 2062, 2027, 2060, 2061, 2070, 2032, 7154, 7153, 11779, 5932, 2976, 4613, 794, 1043, 4632, 4631, 4633, 4404, 2059, 12095, 2040, 2915, 4372, 2069, 11229, 7848, 15383, 10125, 2031, 10875, 12889, 12888, 321, 322, 9551, 9550, 3954, 2020, 2068, 4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552])
    #Driver([6096, 1971, 1972, 768, 952, 11674, 13255, 13253, 13246, 12952, 12953, 13247, 12951, 2025])
    #Driver([9548, 2149, 2197, 9510, 6131, 5908, 803, 14342, 12995, 2064, 16024, 12992, 14332, 13202])
    #Driver([2056, 2055, 2922, 9506, 11344, 766, 4688, 6869, 6872, 3554, 2057, 2058, 8041, 9121, 9546, 7252, 7060, 9553, 5930, 5931, 5929, 2079, 5905, 9527, 4689, 3947, 1563, 9507, 4613, 794, 11775, 11271, 3951, 2062, 2027, 2060, 2061, 2070, 2032, 7154, 7153, 11779, 5932, 2976, 4613, 794, 1043, 4632, 4631, 4633, 4404, 2059, 12095, 2040, 2915, 4372, 2069, 11229, 7848, 15383, 10125, 2031, 10875, 12889, 12888, 321, 322, 9551, 9550, 3954, 2020, 2068, 4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552])
    #Driver([4404, 2059, 12095, 2040, 2915, 4372, 2069, 11229, 7848, 15383, 10125, 2031, 10875, 12889, 12888, 321, 322, 9551, 9550, 3954, 2020, 2068, 4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552])
    #Parallelization Testing
    #Driver([6096])
    #Driver([1971])
    #Driver([6096, 1971, 1972, 768])
    #Driver([6096, 1971])
    #Driver([1972])
    #Parallelization
    #[4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552]
    #Driver([4742, 2039, 3150, 2030, 4743, 5197, 11784, 9552])
    #ObsID_L=ObsID_From_CSV_Query.Read_ObsIDs(Remove_Unarchived=True)
    ObsID_L=ObsID_From_CSV_Query.Read_ObsIDs()
    Driver(ObsID_L)

Main()
