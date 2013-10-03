function photobin,wavl,flux,filterX,filterY,showplot=showplot
;; Bin a spectrum from a given photometry responsivity curve
;; wavl - wavelenght for spectrum
;; flux - flux for spectrum
;; filterX/Y - the filter curve
;; showplot shows the binning process

;; Make sure the filter ends are zero so when extrapolating all fluxes
;; outside the range are weighted as zero
dlam = filterX[1] - filterX[0]
nfilter = n_elements(filterX)
filterX = [filterX[0] - 2E * dlam,filterX[0] - dlam,$
           filterX,filterX[nfilter-1l] + dlam,filterX[nfilter-1l] + 2E * dlam]
filterY = [0E,0E,filterY,0E,0E]

;; interpolate the filter curve
filterInterp = interpol(filterY,filterX,wavl)

plot,filterInterp * flux
binnedValue = total(filterInterp * flux) / total(filterInterp)

if keyword_set(showplot) then begin
   plot,wavl,flux,ystyle=16
   peakval = max(filterInterp,peakInd)
   oplot,[wavl[peakInd]],[binnedValue],psym=2,symsize=4,$
         color=mycol('red')
   
   prevXrange = !x.crange
   plot,wavl,filterInterp,color=mycol('yellow'),$
        /noerase,xrange=prevXrange,$
        xstyle=4+1,ystyle=16+4
   stop
endif

return,binnedValue
end
