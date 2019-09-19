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
pro idl_to_csv

  restore, 'BEHR_results_allsources.sav'
  headers = ['obsid','source','hardmean','softmean','harderrorpos','harderrorneg','softerrorpos','softerrorneg']
  HR1Error=transpose(HR1Error)
  HR2Error=transpose(HR2Error)
  hardErrorPos = hr2error[*,0]
  hardErrorNeg = hr2error[*,1]
  softErrorPos = hr1error[*,0]
  softErrorNeg = hr1error[*,1]
  print, n_elements(hr2mean)
  print, n_elements(hr1mean)
  print, n_elements(hardErrorPos)
  print, n_elements(hardErrorNeg)
  print, n_elements(softErrorPos)
  print, n_elements(softErrorNeg)

  write_csv,'behr_results_allsources.csv',obsid,source,hr2mean,hr1mean,harderrorpos,harderrorneg,softerrorpos,softerrorneg,header=headers
  stop
end
