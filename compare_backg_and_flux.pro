pro compare_backg_and_flux,showplot=showplot,specialInd=specialInd
;; Comapares the flux ratio to the background ratio
;; showplot - shows the correlation
  
  ;; get the compiled spectroscopic data
  restore,'data/specdata.sav'

  ;; get the time info
  restore,'data/timedata.sav'

  ;; divide all time series by the transit model
  ymodel = quadlc(tplot,planetdat.p,planetdat.b_impact,$
                  u1parm,u2parm,planetdat.a_o_rstar)
  divsize = size(DivSpec)
  replicatedmodel = rebin(ymodel,divsize[3],divsize[2],divsize[1])
  rebinmodel = transpose(replicatedmodel,[2,1,0])
  divbycurve = DivSpec / rebinmodel

  nwavs = n_elements(lamgrid)
  slopeArray = fltarr(nwavs)
  offsetArray = fltarr(nwavs)
  back2source = fltarr(nwavs)

  !p.multi = [0,0,2]

  for i=0l,nwavs-1l do begin
     xnorm = backdiv[i,0,*] ;* median(backgrid[i,0,*]) /median(flgrid[i,0,*])
     xnorm = xnorm ;/ median(xnorm)
     back2source[i] = median(backgrid[i,0,*]) /median(flgrid[i,0,*])
     
     ynorm = divbycurve[i,0,*]
     ynorm = ynorm/median(ynorm)
     goodp = where(finite(xnorm) and finite(ynorm))
     if n_elements(goodp) GT 20 then begin
        rfit = robust_linefit(xnorm[goodp],ynorm[goodp],yfit)
        
        slopeArray[i] = rfit[1]
        offsetArray[i] = rfit[0]
     endif
;     specialInd =200
;     specialInd =448
     if n_elements(specialInd) EQ 0 then specialInd =150
     if keyword_set(showplot) and i GE specialInd then begin
        if i EQ specialInd then begin
           !p.multi = [0,0,3]
 
           plot,xnorm,ynorm,$
                xrange=threshold(xnorm,mult=0.3,low=0.02,high=0.98),$
;                xrange=[0.1,1.0],$
                yrange=threshold(ynorm,mult=0.6,low=0.02,high=0.98),psym=4,$
                xtitle='Background ratio',$
                ytitle='Normalzed flux ratio',charsize=2
           xshow = findgen(64)/64E * (!x.crange[1] - !x.crange[0]) + !x.crange[0]
           yshow = rfit[0] + rfit[1] * xshow
           oplot,xshow,yshow,color=mycol('red')
        endif else junk=junk;oplot,xnorm,ynorm,psym=4

     endif

     
  endfor

  plot,back2source,slopeArray,yrange=[-1,0.5],$
;  plot,lamgrid,slopeArray,yrange=[-2,1],$
 ;      xtitle='Background to Source Ratio',charsize=2$
       ytitle='Line slope',charsize=2,psym=4,xrange=[0.03,2],/xlog

  ;; Show the slope of the highlighted index
  oplot,[back2source[specialInd]],[slopeArray[specialInd]],$
        psym=4,color=mycol('red'),thick=2

;  slopeMed = median(slopeArray)
;  oplot,[lamgrid[0],lamgrid[nwavs-1]],replicate(slopeMed,2),color=mycol('lblue')
;  print,'Median slope = ',slopeMed

;  plot,lamgrid,offsetArray,yrange=[0,4],$
  plot,back2source,offsetArray,yrange=threshold(offsetarray,low=0.1,high=0.97,mult=0),$
;       xtitle='Wavelength (um)',$
       xtitle='Background to source ratio',$
;       ytitle='Line fit Offset',charsize=2
       ytitle='Line fit Offset',charsize=2,/xlog,$
       xrange=threshold(back2source,high=0.95,low=0.02,mult=0),psym=4,$
       xstyle=1

  ;; Show the special index
  oplot,[back2source[specialInd]],[offsetArray[specialInd]],$
       psym=4,color=mycol('red'),thick=2

  !p.multi = 0

end
