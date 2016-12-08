function make_tplot,utgrid,timeName=timeName,hr=hr
;; Makes a time array with a JD reference
  
  medTime = median(utgrid)
  refTime = double(round(medTime))
  tplot = utgrid - refTime
  timeName = 'JD - '+string(refTime,format='(I8)')
  if keyword_set(hr) then begin
     tplot = tplot * 24E
     timeName = timeName + ' (hr)'
  endif

  return, tplot
end
