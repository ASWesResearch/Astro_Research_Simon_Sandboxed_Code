"""
This code is the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
"Exploring the High Energy Sky: A Survey of X-ray Sources in Galaxies Within 15 Megaparsecs of the Milky Way".
This thesis was submitted to the faculty of Wesleyan University in partial fulfillment of the requirements for the Degree of Bachelor of Arts
with Departmental Honors in Astronomy.

This version of the code is now modifed by Anthony Santini (Wesleyan 2019 Master's alumnus) for use in current research.


"""
import os
#from ciao_contrib.runtool import *
import numpy as np
import glob
from done_email import when_done


def clean_results(obsid):
    #the following makes all obsids into integers and removes anything globbed that isnt a number (and therefore not an obsid)
    for i in range(0,len(obsid)):
        try:
            obsid[i] = int(obsid[i])
        except ValueError:
            obsid[i] = -1
    obsid.remove(-1)
    obsid.remove(-1) #THIS IS A PRETTY POOR SOLUTION, BUT IT WORKS FOR NOW
    #print len(obsid)
    #for item in obsid:
    #    if item == -1:
    #        obsid.remove(item)
    #print len(obsid)
    #above should be 2 shorter than first length, since there are 2 non-obs id files in that folder
    return obsid

def spec(obsid):
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
    for i in range(1,num+1):

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
        print "Source " + str(i) + " of " + str(num) + " for obsid " + str(obsid) + "..."
        os.system("specextract \
        infile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + ".reg)]' \
        bkgfile='chandra_from_csc/" + str(obsid) + "/primary/" + str(name) + "[sky=region(/Volumes/xray/spirals/trace/" + str(obsid) + "/" + str(i) + "_bkg.reg)]' \
        outroot='chandra_from_csc/" + str(obsid) + "/primary/extracted/extracted_" + obsid + "_" + str(i) + "' \
        correctpsf=yes \
        weight=no \
        clobber=yes \
        verbose=1")


if __name__ == '__main__':

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

    for i in range(0,len(obsid)):
        try:
            #print "Attempting specextract..."
            spec(obsid[i])
        except OSError:
            print "No traced sources in obsid " + str(obsid[i])
            print "specextract failed for obsid " + str(obsid[i])

    when_done("specextract")
