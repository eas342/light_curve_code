function find_shifts,inArray,cutEnds=cutEnds,stopAndshow=stopAndshow
;pro find_shifts
;; Finds the wavelength shifts for a series of spectra
;; Assumes that the X direction is wavelength & Y Direction is
;; spectrum #
;; cutEnds - trims the bottom & top 10% of spectrum before cross-correlating
;; stopAndshow - stops the cross-correlation to show
;;               cross-correlations, peaks and the spectrum

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
;  medspec = inarray[*,0]

  ;; Make a lag array
  nLag = 20l
;  nLag = 100l
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
     ;; Discard all extreme points
     if keyword_set(cutEnds) then begin ;; cut out bottom & top 10%
        badp = where(lindgen(nwavs) LT float(nwavs)*0.1E OR $
                     lindgen(nwavs) GT float(nwavs)*0.9E OR $
                     inArray[*,j] EQ 0 OR $
                     medspec EQ 0,$; OR $
;                     (lindgen(nwavs) GE 220 and lindgen(nwavs) LT 235),$
                     numbad,complement=keepp)
     endif else keepp=lindgen(nwavs)
     numgood = nwavs - numbad
     if numgood LE nlag then shiftArr[j] = 0 else begin
        
        crossCor = c_correlate(inArray[keepp,j],medspec[keepp],lagArray)
;    peakVal = max(crossCor,peakP) 
        PolyFit = poly_fit(lagArray,crossCor,2,yfit=polyVals)
        shiftArr[j] = polyFit[1]/(-2E * polyFit[2])
;     if lagArray[peakP] NE 0 then begin

     endelse
;;     outArray[*,j] = shift_interp(inArray[*,j],shiftArr[j])
        
;        outArray[*,j] = shift_interp(inArray[*,j],shiftArr[j])
        outArray[*,j] = shift_interp(inArray[*,j],shiftArr[j])
;        outArray[*,j] = shift_interp(inArray[*,j],-result[2])

     if j GE 575 and keyword_set(stopAndshow) then begin
        !p.multi = [0,1,2]
;     endif
;     if j GE 0 then begin
        plot,lagArray,crossCor,ystyle=16,$
             xtitle='Shift (px)',ytitle='Cross Cor',psym=2
        oplot,lagArray,PolyVals,color=mycol('green')
        oplot,[shiftArr[j],shiftArr[j]],!y.crange,color=mycol('yellow')
        fitExpr = 'P[0] * SinC(P[1] * (X -P[2])) + P[3] + P[4] *  EXP(-0.5E * ((X - P[2])/P[5])^2)'
        startParams = [1,0.3,0,8,0,3]
        result=mpfitexpr(fitExpr,lagarray,crosscor,fltarr(nlag)+0.1,$
                       startparams,/quiet)
        oplot,lagarray,expression_eval(fitExpr,lagarray,result),color=mycol('blue')
        oplot,result[[2,2]],!y.crange,color=mycol('red')
        legend,['Cross Correlation','Best-Fit Polynomial','Best fit sinc + Gaussian',$
                'Polynomial Peak','Gaussian + Sinc Peak'],$
               linestyle=[0,0,0,0,0],psym=[1,0,0,0,0],color=mycol(['white','green','blue','yellow','red'])

        outArray[*,j] = shift_interp(inArray[*,j],result[2])

        x = findgen(nwavs)
        plot,x,inArray[*,j],ystyle=16,$
             xtitle='Pixel',ytitle='Flux',xrange=[100,400]
        oplot,x,medspec*0.97,color=mycol('green')
        oplot,x,outArray[*,j],color=mycol('red')
        oplot,x,shift_interp(inArray[*,j],shiftArr[j]),color=mycol('yellow')
        legend,['Original','Median Spec','W/ Polynomial Peak','W/ Sinc Peak'],$
               color=mycol(['white','green','red','yellow']),$
               linestyle=[0,0,0,0],/right

        stop
        !p.multi = 0
     endif

  endfor


  return,outArray
end
