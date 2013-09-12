pro analyze_resids,psplot=psplot
  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/power_spectrum/acf_plot'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
  endif


;; Looks at the residuals of the model fit to analyze red noise

  cd,c=currentd
  fileopt = file_search(currentd+'/data/cleaned_tim_ser/*.txt')
  totfiles = n_elements(fileopt)
;  for i=0l,n_elements(fileopt)-1l do begin
  readcol,fileopt[0],phase,fl,flerr,modelfl,resid,format='(F)',skipline=1,/silent
;     readcol,'data/cleaned_tim_ser/timeser_1.43um_.txt',$
;  readcol,'data/cleaned_tim_ser/timeser_0.91um_.txt',$
;          phase,fl,flerr,modelfl,resid

  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
          skipline=1
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor

  t = phase * planetdat.period * 24D * 60D ;; min

  ;; Find autocorrelation function
  np = n_elements(phase)
  steparray = lindgen(np)
  autoC = a_correlate(resid,steparray)
  if keyword_set(min) then begin
     autoX = steparray * (t[1] - t[0])
     autoXtitle = 'Delay (min)'
     autoXrange = [0.4,200]
  endif else begin
     autoX = steparray
     autoXtitle = 'Lag (steps)'
     autoXrange = [1,n_elements(steparray)-1l]
  endelse
  plot,autoX,autoC,$
       xtitle=autoXtitle,$
       ytitle='Autocorrelation',$
       title='1.43um Time Series',$
       xrange=autoXrange

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
     plotprenm = 'plots/power_spectrum/psd_plot'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
  endif

  scargle,t,fl,om,px
  plot,om/(2D * !DPI),px,$
;       xtitle='Frequency (1/min)',$
       xtitle='Frequency (1/min)',/xlog,$
;       ytitle='Power Spectral Density',yrange=[0,5],$
       ytitle='Power Spectral Density',/ylog,$
       title='1.43um Time Series'


  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
     plotprenm = 'plots/power_spectrum/residual_series'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
  endif

  plot,t,resid,xtitle='Time from transit center (s)',$
       ytitle='Flux Residuals (%)',$
       title='1.43um Time Series'
;  oploterr,t,resid,flerr*100E
;  stop
  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif


end
