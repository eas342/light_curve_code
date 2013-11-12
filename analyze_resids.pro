pro analyze_resids,psplot=psplot,showkern=showkern,fast=fast
;; Looks at the residuals of the model fit to analyze red noise
;; also generates, Autocorrelation and Power Spectral Density
;; showkern -- show the covariance kernel w/ best-fit hyper-parameters
;; fast - skips the time series and power spectrum plots to run faster

  cd,c=currentd
  fileopt = file_search(currentd+'/data/cleaned_tim_ser/*.txt')
  totfiles = n_elements(fileopt)
;  for i=0l,n_elements(fileopt)-1l do begin
;     readcol,'data/cleaned_tim_ser/timeser_1.43um_.txt',$
;  readcol,'data/cleaned_tim_ser/timeser_0.91um_.txt',$
;          phase,fl,flerr,modelfl,resid

  if keyword_set(showkern) then begin
     readcol,'radius_vs_wavelength/fit_data_mcmc/09_9Q_X_D0_vs_wavl.txt',$
             format='(F,F,F,F)',wavl1,wavl1size,theta0,theta0Err
     readcol,'radius_vs_wavelength/fit_data_mcmc/10_9Q_X_D1_vs_wavl.txt',$
             format='(F,F,F,F)',wavl2,wavl2size,theta1,theta1Err
  endif

  for wavInd = 0l,totfiles-1l do begin
     readcol,fileopt[wavInd],phase,fl,flerr,modelfl,resid,format='(F)',skipline=1,/silent
     
     startString = 'cleaned_tim_ser/timeser_'
     wavnSt = strpos(fileopt[wavInd],startString)
     wavnSize = strpos(fileopt[wavInd],'_.txt') - wavnSt - strlen(startString)
     wavname = strmid(fileopt[wavInd],wavnSt+strlen(startString),wavnSize)
     
     ;; set the plot
     if keyword_set(psplot) then begin
        set_plot,'ps'
        !p.font=0
        plotprenm = 'plots/power_spectrum/acf_plot_'+wavname
        device,encapsulated=1, /helvetica,$
               filename=plotprenm+'.eps'
        device,xsize=12, ysize=8,decomposed=1,/color
     endif
     
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
     if keyword_set(showkern) then begin
        autoC = a_correlate(resid * 1E-2,steparray,/cov)
        custYrange=[min(autoC),max(autoC) * 1.5E]
     endif else begin
        autoC = a_correlate(resid,steparray)
        custYrange=[0,0]
     endelse
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
          title=wavname+' Time Series',$
          xrange=autoXrange,$
          yrange=custYrange

     if keyword_set(showkern) then begin
        ;; Show the best-fit kenrel
        ;; get the MCMC hyperpameter fit data
        kernX = phase - phase[0]
        kernY = cov_kernel(kernX,theta0[wavInd],theta1[wavInd])
        oplot,autoX,kernY,color=mycol('blue')
        kernY2 = cov_kernel(kernX,theta0[wavInd] - theta0Err[wavInd],$
                            theta1[wavInd])
        oplot,autoX,kernY2,color=mycol('blue'),linestyle=2
        kernY3 = cov_kernel(kernX,theta0[wavInd] + theta0Err[wavInd],$
                            theta1[wavInd])
        oplot,autoX,kernY3,color=mycol('blue'),linestyle=2
;        kernY4 = cov_kernel(kernX,theta0[wavInd],$
;                            theta1[wavInd] + theta1Err[wavInd])
;        oplot,autoX,kernY4,color=mycol('red'),linestyle=2
;        kernY5 = cov_kernel(kernX,theta0[wavInd],$
;                            theta1[wavInd] - theta1Err[wavInd])
;        oplot,autoX,kernY5,color=mycol('red'),linestyle=2

     endif

     if not keyword_set(fast) then begin
        if keyword_set(psplot) then begin
           device, /close
           cgPS2PDF,plotprenm+'.eps'
           spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
           plotprenm = 'plots/power_spectrum/psd_plot_'+wavname
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
             title=wavname+' Time Series'
        
        
        if keyword_set(psplot) then begin
           device, /close
           cgPS2PDF,plotprenm+'.eps'
           spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
           plotprenm = 'plots/power_spectrum/residual_series_'+wavname
           device,encapsulated=1, /helvetica,$
                  filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
        endif
     
        plot,t,resid,xtitle='Time from transit center (s)',$
             ytitle='Flux Residuals (%)',$
             title=wavname+' Time Series'

     endif 

     if keyword_set(psplot) then begin
        device, /close
        cgPS2PDF,plotprenm+'.eps'
        spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
        device,decomposed=0
        set_plot,'x'
        !p.font=-1
     endif
  endfor
  
end
