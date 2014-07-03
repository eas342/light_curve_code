pro plot_amplitudes,psplot=psplot,secondary=secondary
;; Plots the amplitudes of the fitted profiles
;; this is an alternative to IRAF extraction
;; psplot -- make a postscript plot
;; secondary -- use secondary eclipse timing (centered around 0.5
;;              instead of 0.0)

  restore,'data/used_date.sav'

  restore,'data/state_parameters/widths/widths_'+$
          specfileListNamePrefix+'.sav'

  plot_tim_ser,secondary=secondary
  restore,'data/timedata.sav'
  
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/spec_t_series/amplitudes'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     device,xsize=15, ysize=10,decomposed=1,/color
     !p.thick=3.5
     !x.thick=2
     !y.thick=2
  endif

  goodp = where(staramps[*,1] NE 0)
  ntime = n_elements(staramps[*,0])
  y = fltarr(ntime) * !values.f_nan
  y[goodp] = staramps[goodp,0] / staramps[goodp,1]
  y = y/median(y)
  plot,tplot,y,psym=4,yrange=threshold(y),$
       xtitle='Orbital Phase',$
       ytitle='Normalized Amplitude',$
       xstyle=1

  if keyword_set(psplot) then begin
     device, /close
;     cgPS2PDF,plotprenm+'.eps'
;     spawn,'convert -density 300% '+plotprenm+'.pdf'+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick=1
     !x.thick=1
     !y.thick=1
  endif
end
