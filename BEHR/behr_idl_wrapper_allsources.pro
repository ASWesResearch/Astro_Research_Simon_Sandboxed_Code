;"""
;This code is the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
;"Exploring the High Energy Sky: A Survey of X-ray Sources in Galaxies Within 15 Megaparsecs of the Milky Way".
;This thesis was submitted to the faculty of Wesleyan University in partial fulfillment of the requirements for the Degree of Bachelor of Arts
;with Departmental Honors in Astronomy.
;
;This version of the code is now modifed by Anthony Santini (Wesleyan 2019 Master's alumnus) for use in current research.
;
;
;"""
pro behr_idl_wrapper

  ;This is currently working for not_csc data included in the counts_info.csv document. The plots must be done afterwards, for some reason.
  ;Must be modified to include data that will come from the not_csc observations. Must also include checking against list of good not_csc sources.


  lun=1
  ;openr, lun,'counts_info_all_clean.csv'
  openr, lun,'counts_info.csv'
  line=''
  Obsid=[]
  Source=[]
  SoftC=[]
  MedC=[]
  HardC=[]
  SoftBKG=[]
  MedBKG=[]
  HardBKG=[]

  readf, lun, line

  while ~ eof(lun) do begin
    readf,lun, line
    l = strsplit(line, ',', /extract)
    Obsid=[Obsid,l[0]]
    Source=[Source,l[1]]
    SoftC=[SoftC,double(l[2])]
    MedC=[MedC,double(l[3])]
    HardC=[HardC,double(l[4])]
    SoftBKG=[SoftBKG,double(l[5])]
    MedBKG=[MedBKG,double(l[6])]
    HardBKG=[HardBKG,double(l[7])]
  endwhile
  close, lun

;  lun=1
;  openr, lun,'goodsources_nodupes.txt'
;  line=''
;  goodObsid=[]
;  goodSource=[]
;  readf, lun, line
;  while ~ eof(lun) do begin
;    readf,lun, line
;    l = strsplit(line, /extract)
    ;print, l
;    goodObsid=[goodObsid,l[0]]
;    goodSource=[goodSource,l[1]]
;  endwhile
;  close, lun

  ;print, source[0]


  HR1Mean=[]
  HR1Upper=[]
  HR1Lower=[]
  HR1Error=[]
  HR2Mean=[]
  HR2Upper=[]
  HR2Lower=[]
  HR2Error=[]

  for i=0,n_elements(SoftC)-1 do begin
    print, Obsid[i], " ", Source[i]

    softCurr = round(SoftC[i])
    medCurr = round(MedC[i])
    hardCurr = round(HardC[i])
    softbkgCurr = round(SoftBKG[i])
    medbkgCurr = round(MedBKG[i])
    hardbkgCurr = round(HardBKG[i])

    print, i,softCurr,medCurr,hardCurr,softbkgCurr,medbkgCurr,hardbkgCurr

    behr3=colbehr(softCurr,medCurr,hardCurr,Sbkg=softbkgCurr,Mbkg=medbkgCurr,Hbkg=hardbkgCurr,Sarea=3,Marea=3,Harea=3)

    HR1Mean=[HR1Mean,-behr3.hr1.mean]
    HR1Upper=[HR1Upper,abs(behr3.hr1.upperbound-behr3.hr1.mean)]
    HR1Lower=[HR1Lower,abs(behr3.hr1.lowerbound-behr3.hr1.mean)]

    HR2Mean=[HR2Mean,-behr3.hr2.mean]
    HR2Upper=[HR2Upper,abs(behr3.hr2.upperbound-behr3.hr2.mean)]
    HR2Lower=[HR2Lower,abs(behr3.hr2.lowerbound-behr3.hr2.mean)]
  endfor

  HR1Error=[[HR1Lower],[HR1Upper]]
  HR2Error=[[HR2Lower],[HR2Upper]]
  HR1Error=transpose(HR1Error)
  HR2Error=transpose(HR2Error)

  save, /variables, filename='BEHR_results_allsources.sav'

;  graphic = ERRORPLOT(hr2mean,hr1mean,hr2error,hr1error,LINESTYLE=' ',TITLE='Hardness Ratios')
;  ax = graphic.AXES
;  ax[0].TITLE='(H-M)/(H+M)'
;  ax[1].TITLE='(M-S)/(M+S)'
;  graphic2 = PLOT(hr2mean, hr1mean,linestyle=' ',symbol='o',color='g',/OVERPLOT)

end
