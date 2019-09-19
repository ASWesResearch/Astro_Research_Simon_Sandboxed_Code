"""
This code is the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
"Exploring the High Energy Sky: A Survey of X-ray Sources in Galaxies Within 15 Megaparsecs of the Milky Way".
This thesis was submitted to the faculty of Wesleyan University in partial fulfillment of the requirements for the Degree of Bachelor of Arts
with Departmental Honors in Astronomy.

This version of the code is now modifed by Anthony Santini (Wesleyan 2019 Master's alumnus) for use in current research.


"""
# Extracts counts from the from_csc sources processed by the raytracer and corrects for area. Creates counts_info.csv

def clean_results(obsid):
    #the following makes all obsids into integers and removes anything globbed that isnt a number (and therefore not an obsid)
    for i in range(0,len(obsid)):
        try:
            obsid[i] = int(obsid[i])
        except ValueError:
            obsid[i] = -1
    obsid.remove(-1)
    obsid.remove(-1) #THIS IS A PRETTY POOR SOLUTION, BUT IT WORKS FOR NOW
    return obsid

def spec(obsid):
    #Not used in this program -- can probably be deleted
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
    name = glob.glob("chandra_from_csc/" + obsid + "/primary/*_evt2.fits")
    print name
    print "Shortening name..."
    name = name[0][-24:]
    print name
    print "Getting system path..."
    path = '/Volumes/xray/spirals/trace/' + obsid + '/'
    print "Getting number of files in the current obsid's directory..."
    num_files = len([f for f in os.listdir(path)
                if os.path.isfile(os.path.join(path, f))])
    print str(num_files) + " files found..."
    num = int((num_files-4)/5)
    print str(num) + " sources found..."

    #######################################################
    ###### Get the Cycle and Instrument for the Data ######
    #######################################################

    print "Getting cycle..."
    date = sp.check_output("dmkeypar /Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + name + " DATE-OBS echo+", shell=True)
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

    for i in range(1, num+1):
        #This was taken from filter_background_sources -- will use it to get the pointing which will tell us about the instrument
        #Finds x and y for a specific source -- then we can determine which chip it is on, and therefore what "instrument" to use

        print "\nReading in for file number " + str(i) + "..."
        f = open("/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + ".reg")
        file_read = f.read()
        first = file_read.split("(")
        first = first[1]
        first = first.split(",")
        x = float(first[0])
        y = float(first[1])
        f.close()
        print "Got x and y values..."

        print "punlearn dmcoords..."
        os.system("punlearn dmcoords")
        print "running dmcoords..."
        #print "dmcoords '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "' option=sky x=" + str(x) + " y=" + str(y) + " celfmt=deg"
        os.system("dmcoords '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "' option=sky x=" + str(x) + " y=" + str(y) + " celfmt=deg")
        print "Getting Chip ID..."
        chip = sp.check_output("pget dmcoords chip_id", shell=True)
        chip = int(chip)
        print "Got Chip ID: " + str(chip)

        if chip == 5 or chip == 7:
            inst = 's'
        else:
            inst = 'i'

        print inst

        #########################
        ###### Soft X-rays ######
        #########################

        print "Extracting soft X-rays for source " + str(i) + " in obsid " + obsid + "..."

        os.system("dmextract \
        infile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=300:1000][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + ".reg)]\" \
        outfile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_softcounts.fits\" \
        opt=generic \
        bkg=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=300:1000][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + "_bkg.reg)]\" \
        clobber=no")

        #Getting counts
        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_softcounts.fits[cols counts]'")
        count = sp.check_output("pget dmstat out_sum", shell=True)
        count = effectiveAreaCorrect(count, cycle, 0, inst)
        rawCountsSoft.append(count)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_softcounts.fits[cols area]'")
        area = sp.check_output("pget dmstat out_sum", shell=True)
        areaSoft.append(area)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_softcounts.fits[cols bg_counts]'")
        bgCount = sp.check_output("pget dmstat out_sum", shell=True)
        bgCount = effectiveAreaCorrect(bgCount, cycle, 0, inst)
        bkgCountsSoft.append(bgCount)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_softcounts.fits[cols bg_area]'")
        bgArea = sp.check_output("pget dmstat out_sum", shell=True)
        bkgAreaSoft.append(bgArea)


        ########################
        ###### Med X-rays ######
        ########################

        print "Extracting medium X-rays for source " + str(i) + " in obsid " + obsid + "..."

        os.system("dmextract \
        infile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=1000:2100][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + ".reg)]\" \
        outfile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_medcounts.fits\" \
        opt=generic \
        bkg=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=1000:2100][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + "_bkg.reg)]\" \
        clobber=no")

        #Getting counts
        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_medcounts.fits[cols counts]'")
        count = sp.check_output("pget dmstat out_sum", shell=True)
        count = effectiveAreaCorrect(count, cycle, 1, inst)
        rawCountsMed.append(count)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_medcounts.fits[cols area]'")
        area = sp.check_output("pget dmstat out_sum", shell=True)
        areaMed.append(area)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_medcounts.fits[cols bg_counts]'")
        bgCount = sp.check_output("pget dmstat out_sum", shell=True)
        bgCount = effectiveAreaCorrect(bgCount, cycle, 1, inst)
        bkgCountsMed.append(bgCount)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_medcounts.fits[cols bg_area]'")
        bgArea = sp.check_output("pget dmstat out_sum", shell=True)
        bkgAreaMed.append(bgArea)


        #########################
        ###### Hard X-rays ######
        #########################

        print "Extracting hard X-rays for source " + str(i) + " in obsid " + obsid + "..."

        os.system("dmextract \
        infile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=2100:7500][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + ".reg)]\" \
        outfile=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_hardcounts.fits\" \
        opt=generic \
        bkg=\"/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/" + str(name) + "[energy=2100:7500][sky=region(/Volumes/xray/spirals/trace/" + obsid + "/" + str(i) + "_bkg.reg)]\" \
        clobber=no")

        #Getting counts
        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_hardcounts.fits[cols counts]'")
        count = sp.check_output("pget dmstat out_sum", shell=True)
        count = effectiveAreaCorrect(count, cycle, 2, inst)
        rawCountsHard.append(count)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_hardcounts.fits[cols area]'")
        area = sp.check_output("pget dmstat out_sum", shell=True)
        areaHard.append(area)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_hardcounts.fits[cols bg_counts]'")
        bgCount = sp.check_output("pget dmstat out_sum", shell=True)
        bgCount = effectiveAreaCorrect(bgCount, cycle, 2, inst)
        bkgCountsHard.append(bgCount)

        os.system("dmstat '/Volumes/xray/simon/chandra_from_csc/" + obsid + "/primary/extracted_counts_info/" + obsid + "_" + str(i) + "_hardcounts.fits[cols bg_area]'")
        bgArea = sp.check_output("pget dmstat out_sum", shell=True)
        bkgAreaHard.append(bgArea)

        obs.append(obsid)
        src.append(i)

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



if __name__ == '__main__':

    import os
    import numpy as np
    import glob
    import subprocess as sp
    from done_email import when_done
    from ciao_contrib.runtool import *
    import csv

    print "Reading in obsids... "
    obsid_path = glob.glob("chandra_from_csc/*")

    obsid = []
    for i in range(0,len(obsid_path)):
        #print obsid_path[i]
        obsid.append(obsid_path[i][17:])
    print "...done!"


    obsid = clean_results(obsid)


    ##############################
    ###### Start Processing ######
    ##############################


    print "All obsid names processed."

    obsidData = ["OBSID"]
    sourceData = ["SOURCE"]
    rawCountsSoft = ["RAW_COUNTS[.3-1]"]
    rawCountsMed = ["RAW_COUNTS[1-2.1]"]
    rawCountsHard = ["RAW_COUNTS[2.1-7.5]"]
    areaSoft = ["AREA[.3-1]"]
    areaMed = ["AREA[1-2.1]"]
    areaHard = ["AREA[2.1-7.5]"]
    bkgCountsSoft = ["BKG_COUNTS[.3-1]"]
    bkgCountsMed = ["BKG_COUNTS[1-2.1]"]
    bkgCountsHard = ["BKG_COUNTS[2.1-7.5]"]
    bkgAreaSoft = ["BKG_AREA[.3-1]"]
    bkgAreaMed = ["BKG_AREA[1-2.1]"]
    bkgAreaHard = ["BKG_AREA[2.1-7.5]"]


    for i in range(0,len(obsid)):
        try:
            print "\nAttempting flux retrieval..."
            obs, src, countS, countM, countH, areaS, areaM, areaH, bgcountS, bgcountM, bgcountH, bgareaS, bgareaM, bgareaH = src_flux(obsid[i])

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
                areaSoft.append(areaS[j][:-1])
                areaMed.append(areaM[j][:-1])
                areaHard.append(areaH[j][:-1])
                bkgCountsSoft.append(bgcountS[j])
                bkgCountsMed.append(bgcountM[j])
                bkgCountsHard.append(bgcountH[j])
                bkgAreaSoft.append(bgareaS[j][:-1])
                bkgAreaMed.append(bgareaM[j][:-1])
                bkgAreaHard.append(bgareaH[j][:-1])


        except OSError:
            #print "No traced sources in obsid " + str(obsid[i])
            print "\nFlux extraction failed for obsid " + str(obsid[i])

    zippedList = zip(obsidData, sourceData, rawCountsSoft, rawCountsMed, rawCountsHard, areaSoft, areaMed, areaHard, bkgCountsSoft, bkgCountsMed, bkgCountsHard, bkgAreaSoft, bkgAreaMed, bkgAreaHard)


    # print len(obsidData)
#
    # print obsidData
    # for i in range(0,len(obsidData)):
        # print zippedList[i]

    with open('counts_info.csv', 'wb') as csvfile:
        cWriter = csv.writer(csvfile, delimiter=',')
        for i in range(0,len(obsidData)):
            cWriter.writerow(zippedList[i])


    #when_done("specextract")