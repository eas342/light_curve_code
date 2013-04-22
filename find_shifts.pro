function find_shifts,inArray
;pro find_shifts
;; Finds the wavelength shifts for a series of spectra
;; Assumes that the X direction is wavelength & Y Direction is
;; spectrum #
  restore,'data/specdata.sav'

;  individual = 1
;  inArray = transpose(flgrid[*,individual-1,*],[0,2,1])  
  
  dims = size(inArray)
  nwavs = dims[1]
  nspec = dims[2]

  medspec = fltarr(nwavs)
  ;; Get a median spectrum
  for i=0l,nwavs-1l do begin
     medspec[i] = median(inArray[i,*])
  endfor

  ;; Make a lag array
  nLag = 20l
  lagArray = lindgen(nLag) - nLag/2l

  ;; Zero out all NaNs
  badp = where(finite(medspec) EQ 0)
  medspec[badp] = 0.0E
  badp = where(finite(inArray) EQ 0)
  inArray[badp] = 0.0E

  ;; Shift Array & New image
  shiftArr = fltarr(nspec)
  outArray = fltarr(nwavs,nspec)

  ;; Cross-Correlate to find the shifts
  for j=0l,nspec-1l do begin
     crossCor = c_correlate(inArray[*,j],medspec,lagArray)
;    peakVal = max(crossCor,peakP) 
     PolyFit = poly_fit(lagArray,crossCor,2,yfit=polyVals)
     shiftArr[j] = polyFit[1]/(-2E * polyFit[2])
;     if lagArray[peakP] NE 0 then begin
;        plot,lagArray,crossCor,ystyle=16
;        oplot,lagArray,PolyVals,color=mycol('green')
;        oplot,[shiftArr[j],shiftArr[j]],!y.crange,color=mycol('yellow')
;     endif
     outArray[*,j] = shift_interp(inArray[*,j],shiftArr[j])

  endfor


  return,outArray
end
