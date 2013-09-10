pro prepare_filter_curve
;; Reads in a filter curve and prepares for display & analysis

  readcol,'../calculations/zprime_transmission/zprime_response.txt.csv',skipline=1,$
          wavel,trans

  plot,wavel,trans,$
       xtitle='Wavelength (um)',$
       ytitle='Transmission (%)'

  ;; Find the 1st moment
  firstmom = total(wavel * trans)/total(trans)
  secondmom = sqrt(total((wavel - firstmom)^2 * trans)/total(trans))

  print,'First Moment = ',firstmom,' um'
  print,'sqrt(Second Moment) = ',secondmom

  oplot,[firstmom,firstmom],!y.crange,color=mycol('yellow')
  oplot,[firstmom,firstmom] + secondmom,!y.crange,color=mycol('yellow'),linestyle=1
  oplot,[firstmom,firstmom] - secondmom,!y.crange,color=mycol('yellow'),linestyle=1

  startWav = 0.834E ;; um Half max, left
  endwav = 0.910E ;; um Half max, right
  midwav = (startWav + endwav)/2E
  FWHM = endwav - startWav

  print,'Middle of half widths (um)= ',midwav
  print,'FWHM (um)= ',FWHM

  oplot,[startWav,startWav],!y.crange,color=mycol('lblue')
  oplot,[endWav,endWav],!y.crange,color=mycol('lblue')

end
