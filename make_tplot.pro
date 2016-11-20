function make_tplot,utgrid,timeName=timeName
;; Makes a time array with a JD reference
  
  medTime = median(utgrid)
  refTime = double(round(medTime))
  tplot = utgrid - refTime
  timeName = 'JD - '+string(refTime,format='(I8)')

  return, tplot
end
