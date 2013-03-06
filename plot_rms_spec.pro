pro plot_rms_spec,psplot=psplot,tryclean=tryclean,saveclean=saveclean,$
                  removelinear=removelinear
;; Plots the RMS along the time series for each wavelength in the
;; spectrum
;; psplot -- makes a postscript plot of the RMS spectrum
;; tryclean - tries to clean up the spectrum on a wavelength by
;;            wavelength basis
;; saveclean - saves the white light curve for analysis by plot_tim_ser
;; removelinear -- fit each curve to a line first before finding the RMS

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/rms_spectrum'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif


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
  ntime = n_elements(tplot)
  sigarray = fltarr(nwavs)
  ;; go through each wavelength and find robust sigma in percentage
  cleanedcurve = fltarr(nwavs,ntime)
  cleanfactor = 4E
  for i=0l,nwavs-1l do begin
     sigarray[i] = robust_sigma(divbycurve[i,0,*])/median(divbycurve[i,0,*])
  endfor

  if keyword_set(tryclean) or keyword_set(removelinear) then begin
     for i=0l,nwavs-1l do begin
        if keyword_set(removelinear) then cleanedcurve[i,*] = divbycurve[i,0,*] else begin
           cleanedcurve[i,*] = DivSpec[i,0,*]
        endelse
        badp = where(abs(divbycurve[i,0,*] - median(divbycurve[i,0,*])) GT sigarray[i] * divbycurve[i,0,*] * 6E)
        if badp NE [-1] then cleanedcurve[i,badp] = !values.f_nan
     endfor
  endif

  if keyword_set(removelinear) then begin
     for i=0l,nwavs-1l do begin
        ;; divide by the line to flatten out
        goodp = where(finite(cleanedcurve[i,*]) EQ 1,ngoodp)
        if ngoodp GT 10 then begin
           rlinefit = robust_linefit(tplot[goodp],cleanedcurve[i,goodp],yfit)
           yflat = cleanedcurve[i,goodp] / yfit
           ;; You don't have to divide by the median because
           ;; the division by the line fit already does that
           sigarray[i] = robust_sigma(yflat)
        endif else sigarray[i] = !values.f_nan
        
     endfor
  endif

  plot,lamgrid,sigarray*100E,$
       xtitle='Wavelength (um)',$
       ytitle='Fractional Error (%)',yrange=[0,10]
  oplot,lamgrid,divspecE[*,0,0]/divspec[*,0,0]*100E,color=mycol('blue'),$
        linestyle=2

  legend,['Robust Sigma','Photon Error'],color=mycol(['black','blue']),$
         linestyle=[0,2]

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

  if keyword_set(tryclean) then begin
;     wait,1
     
     ;; Make a special broadband light curve
     specialpt = where(sigarray LT 0.01,nspecial)

     ;; Weight the points by the rms
;     weights=1/sigarray^2
     weights=1/sigarray^2
     ;; discount wavelengths and/or pixels that vary wildly
;     badp = where(sigarray GT 0.01 OR lamgrid LT 0.87 OR lamgrid GT 2.4)
;     badp = where(sigarray GT 0.20)
     if badp NE [-1] then weights[badp] = 0.0E
     weightsCopy = rebin(weights,nwavs,ntime)
     ;; Set all non-finite values of the cleaned curve to have NANS in the weights
     nonfinite = where(finite(cleanedcurve) EQ 0)
     weightsCopy[nonfinite] = !values.f_nan
     ;; weighted sum normalized by the weights for that spectrum
     ;; this accounts for missing data (that has NANs)
     combinedpt = total(cleanedcurve * weightsCopy,1,/nan)/total(weightsCopy,1,/nan)
     
     

     combinedpt2 = total(cleanedcurve[specialpt,*],1)/float(nspecial)
     plot,tplot,combinedpt2,ystyle=16,psym=4
;     plot,tplot,combinedpt2,ystyle=16,psym=4,yrange=[0.40,0.45]
;     oplot,tplot,combinedpt+0.01,psym=5,color=mycol('green')
;     plot,tplot,combinedpt,psym=5,color=mycol('green'),ystyle=16


     ;try out straight avg over wavelength range
     goodrange = where(lamgrid GE 0.9 and lamgrid LT 2.4,ngoodrange)
     combinedpt2 = total(cleanedcurve[goodrange,*],1)/float(ngoodrange)

;     plot,tplot,combinedpt2,ystyle=16,psym=4
;
     if keyword_set(saveclean) then begin
        
     endif
  endif



end
