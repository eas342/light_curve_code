pro molec_wav_regions,psplot=psplot
;psplot -- saves a postscript plot

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'molec_inputs/acetylene'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif

  restore,'molec_inputs/acetylene_absorb_spec500.sav'
  wavl = 1E/wavn3 * 1E4
  Y = absAcetylt500
  molecName = 'C!D2!NH!D2!N'
  molecFileName = 'C2H2'

  plot,wavl,Y,/nodata,xrange=[0.8,2.45],$
       xtitle='Wavelength (um)',ytitle='Absorption'
  smoothY = smooth(Y,2000)
  oplot,wavl,smoothY,color=mycol('blue')

  nwavs = n_elements(wavl)


  threshhold = 0.01
  molecWavs = where(smoothY GT threshhold and wavl LT 2.4)
  AbsArray = lonarr(nwavs)
  AbsArray[molecWavs] = 1l
;  oplot,wavl,AbsArray/2E,color=mycol('red') ;; accepted wavelengths
  ;; Make sure there isn't a bin starting at the last point
  AbsArray[nwavs-1l] = 0l

  shiftedArray = shift(AbsArray,1)
  Dshifted = shiftedArray - AbsArray
;  oplot,wavl,Dshifted

  starts = wavl[where(Dshifted EQ 1)]
  Ends = wavl[where(Dshifted EQ -1)]
  npieces = n_elements(starts)
  assert,npieces,'=',n_elements(Ends),"Starts and ends don't match up"

  binloc = 0.5 ;; what the Y value of the bins should be
  for i=0l,npieces-1l do begin
     oplot,[starts[i],Ends[i]],[binloc,binloc],color=mycol('purple'),thick=15
  endfor
  
  legend,['Smoothed '+molecname+' abs','Bin Wavelengths'],$
         linestyle=[0,0],color=mycol(['blue','purple']),$
         thick=[1,15]

  ;; Save the start and ends of bins
  forprint,starts,Ends,comment='#'+molecFilename+' Wavelengths start (um) end (um)',$
           textout='molec_inputs/'+molecfilename+'_bin_locations.txt'

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif


  ;; get the molecular absorption spectrum

end
