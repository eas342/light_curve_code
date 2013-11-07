pro compile_corot
;; Makes a CoRoT light curve from archival CoRoT data


  restore,'../spacecraft_data/corot_specdata.sav'

  ;; Adjust the older ut dates because the ephemeris older ephemeris
  ;; is slightly off
  utgrid = utgrid + 0.005E

save,bingrid,binfl,binflE,$
     Nwavbins,binsizes,$
     airmass,altitude,$
     utgrid,wavname,$
     filename='data/specdata.sav'
  


end
