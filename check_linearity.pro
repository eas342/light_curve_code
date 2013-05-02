pro check_linearity,psplot=psplot,starN=starN,wavN=wavN
;; checks the time series for detector non-linearity by comparing star
;; ratios to total fluxes
;; psplot -- makes a postscript plot
;; starN -- chooses a star number
;; wavN -- chooses which wavelength bin to look at

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/linearity/star_ratio'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     device,xsize=14, ysize=10,decomposed=1,/color
  endif

  ;; Get flux data
  restore,'data/specdata.sav'
  ;; Get timing data
  restore,'data/timedata.sav'

  if n_elements(starN) EQ 0 then starN = 1 ;default star number

  ;; find the out of transit points
  offp = where(tplot LT hstart OR tplot GT hend,complement=inpt)

  sNames = ['Host Star','Reference Star']

  if n_elements(wavN) EQ 0 then wavN = 0 ;; which wavelength index

  ;; Choose ranges
  robustSig = robust_sigma(binfl[wavN,offp])
  medianVal = median(binfl[wavN,offp])
  myYrange = [medianVal - 10E * RobustSig,medianVal + 4E * RobustSig]
  robustSig = robust_sigma(binind[wavN,starN,offp])
  medianVal = median(binind[wavN,starN,offp])
  myXrange = [medianVal - 4E * RobustSig,medianVal + 4E * RobustSig]
  
  ;; Find out the wavelength bin name
  mytitle = string(bingrid[wavN],format='(F6.3)')+$
            ' to '+string(bingrid[wavN]+binsizes[wavN],format='(F6.3)')+$
            'um Flux Bin'

  plot,binind[wavN,starN,offp],binfl[wavN,offp],$
       xtitle=sNames[starN]+' Flux (e!E-!N)',$
       ytitle='Host Star / Reference Star Ratio',$
       psym=2,xrange=myXrange,yrange=myYrange,$
       title=mytitle
  
  oplot,binind[wavN,starN,inpt],binfl[wavN,inpt],$
        color=mycol('blue'),psym=4

  legend,['Out of Transit','In Transit'],$
         psym=[2,4],color=mycol(['black','blue']),$
         /right

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 200% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

end
