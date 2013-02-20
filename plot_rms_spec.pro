pro plot_rms_spec,psplot=psplot,tryclean=tryclean
;; Plots the RMS along the time series for each wavelength in the
;; spectrum
;; psplot -- makes a postscript plot of the RMS spectrum
;; tryclean - tries to clean up the spectrum on a wavelength by
;;            wavelength basis

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

  if keyword_set(tryclean) then begin
     for i=0l,nwavs-1l do begin
        cleanedcurve[i,*] = DivSpec[i,0,*]
        badp = where(abs(divbycurve[i,0,*] - median(divbycurve[i,0,*])) GT sigarray[i] * divbycurve[i,0,*] * 4E)
        if badp NE [-1] then cleanedcurve[i,badp] = !values.f_nan
     endfor
  endif


  plot,lamgrid,sigarray*100E,$
       xtitle='Wavelength (um)',$
       ytitle='Robust Sigma (%)',yrange=[0,10]

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

  if keyword_set(tryclean) then begin
     wait,1
     
     ;; Make a special broadband light curve
     specialpt = where(sigarray LT 0.01,nspecial)
     combinedpt = total(cleanedcurve[specialpt,*],1)/float(nspecial)
     plot,tplot,combinedpt,ystyle=16,psym=4

  endif



end
