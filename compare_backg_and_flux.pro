pro compare_backg_and_flux,showplot=showplot
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

  for i=0l,nwavs-1l do begin
     xnorm = backdiv[i,0,*] * median(backgrid[i,0,*]) /median(flgrid[i,0,*])
     xnorm = xnorm
     ynorm = divbycurve[i,0,*]
     ynorm = ynorm/median(ynorm)
     goodp = where(finite(xnorm) and finite(ynorm))
     if n_elements(goodp) GT 20 then begin
        rfit = robust_linefit(xnorm[goodp],ynorm[goodp],yfit)
        
        slopeArray[i] = rfit[1]
     endif
     if keyword_set(showplot) and i EQ 70 then begin
        !p.multi = [0,0,2]
        plot,xnorm,ynorm,$
             xrange=median(xnorm) + [-1E,1E] * robust_sigma(xnorm)*3E,$
             yrange=[0.95,1.05],psym=4
        xshow = findgen(64)/64E * (!x.crange[1] - !x.crange[0]) + !x.crange[0]
        yshow = rfit[0] + rfit[1] * xshow
        oplot,xshow,yshow,color=mycol('red')
     endif

     
  endfor

  plot,lamgrid,slopeArray,yrange=[-2,1]

  ;; Show the median slope
  slopeMed = median(slopeArray)
  oplot,[lamgrid[0],lamgrid[nwavs-1]],replicate(slopeMed,2),color=mycol('lblue')
  print,'Median slope = ',slopeMed
  !p.multi = 0

end
