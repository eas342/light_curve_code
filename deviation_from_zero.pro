pro deviation_from_zero
;; Finds out how many sigma from zero the spectrum is (meant for the
;; control night)

  radfile = 'radius_vs_wavelength/radius_vs_wavl.txt'

  readcol,radfile,wavl,wavlsize,rad,rade,skipline=1,format='(F,F,F,F)'

  dev = rad/rade
  print,'dev in sigma= ',dev
end  
  
