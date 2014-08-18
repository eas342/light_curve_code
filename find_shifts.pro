function find_shifts,inArray,cutEnds=cutEnds,stopAndshow=stopAndshow,$
                     masterspec=masterspec
;pro find_shifts
;; Finds the wavelength shifts for a series of spectra
;; Assumes that the X direction is wavelength & Y Direction is
;; spectrum #
;; cutEnds - trims the bottom & top 10% of spectrum before cross-correlating
;; stopAndshow - stops the cross-correlation to show
;;               cross-correlations, peaks and the spectrum
;; masterspec -- an input master spectrum to cross-correlate,
;;               otherwise it uses a median spectrum

  restore,'data/specdata.sav'

;  individual = 1
;  inArray = transpose(flgrid[*,individual-1,*],[0,2,1])  
  
  dims = size(inArray)
  nwavs = dims[1]
  nspec = dims[2]

  ;; Eliminate bad pixels and Nans
  for i=0l,nwavs-1l do begin
     ;; check if whole row is Nans, in that case make it zeros
     finites = finite(inArray[i,*])
     if total(finites) LT nspec - 12l then begin
        inArray[i,*] = 0E
     endif else begin
        rsigma = robust_sigma(inArray[i,*])
        if rsigma EQ [-1] then stop
        medAlong = median(inArray[i,*])
        badp = where(finites EQ 0 OR $
                     abs(inArray[i,*] - medAlong) GT 10E * rsigma)
        if badp NE [-1] then inArray[i,badp] = medAlong
     endelse
  endfor

  ;; Filter the array to get rid of broad features
  finArray = convol(inArray,digital_filter(0.03,0.09,50,25))

  if keyword_set(masterspec) then begin
     badp = where(finite(masterspec) EQ 0)
     cmasterspec = masterspec 
     if badp NE [-1] then cmasterspec[badp] = 0E ;; cleaned masters spec
     medspec = convol(cmasterspec,digital_filter(0.03,0.09,50,25))
  endif else begin
     ;; Get a median spectrum
     medspec = median(finarray,dimension=2)
  endelse

  ;; Make a lag array
  nLag = 26l
  lagArray = lindgen(nLag) - nLag/2l
  fitWidth = 8l ;; size of region used for parabolic fit

  ;; Make a shift Array & New image
  shiftArr = fltarr(nspec)
  outArray = fltarr(nwavs,nspec)

  ;; Cross-Correlate to find the shifts
  for j=0l,nspec-1l do begin
     ;; Discard all extreme points
     if keyword_set(cutEnds) then begin ;; cut out bottom & top 10%
        badp = where(lindgen(nwavs) LT float(nwavs)*0.1E OR $
                     lindgen(nwavs) GT float(nwavs)*0.9E OR $
                     finArray[*,j] EQ 0 OR $
                     medspec EQ 0,$; OR $
;                     (lindgen(nwavs) GE 220 and lindgen(nwavs) LT 235),$
                     numbad,complement=keepp)
     endif else keepp=lindgen(nwavs)
     numgood = nwavs - numbad
     if numgood LE nlag then shiftArr[j] = 0 else begin
        
        crossCor = c_correlate(finArray[keepp,j],medspec[keepp],lagArray)
        peakVal = max(crossCor,peakInd)
        peakP = lagarray[peakInd]

        if abs(peakP) GE nLag/2l - fitWidth/2l then begin
           shiftArr[j] = 0
           print,"Warning, cross correlation peak outside of range"
        endif else begin
           fitPoints = where(lagArray GT peakP-5l and lagarray LT peakP + 5l)
           PolyFit = poly_fit(lagArray[fitpoints],crossCor[fitpoints],2,yfit=polyVals)
           shiftArr[j] = polyFit[1]/(-2E * polyFit[2])
;           plot,lagarray,crosscor
;           oplot,lagarray,PolyFit[0] + PolyFit[1] * lagarray + PolyFit[2] * lagarray^2,color=mycol('yellow')
;           oplot,[shiftArr[j],shiftArr[j]],!y.crange,color=mycol('red')
;           stop

        endelse
     endelse

     outArray[*,j] = shift_interp(inArray[*,j],shiftArr[j])
     if keyword_set(showEach) then begin
        !p.multi=[0,1,2]
        plot,medspec
        oplot,finArray[*,j],color=mycol('red')
        plot,medspec,title='After Shift'
        oplot,shift(finArray[*,j],round(shiftArr[j])),color=mycol('red')
        !p.multi=0
        stop
     endif

     if j GE 350 and keyword_set(stopAndshow) then begin
        !p.multi = [0,1,2]
        plot,lagArray,crossCor,ystyle=16,$
             xtitle='Shift (px)',ytitle='Cross Cor',psym=2
        oplot,lagArray[fitpoints],PolyVals,color=mycol('green')
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

        !p.multi = 0
        stop
     endif

  endfor

  ;; save the shifted array
  save,shiftArr,filename='data/wavelength_shifts/temp_shift_list.sav'

  return,outArray
end
