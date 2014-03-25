pro compare_width_and_flux,fast=fast,allpoints=allpoints,noLine=noLine,$
                           psplot=psplot
;; Compares with spatial profile widths to flux
;; fast -- skips the initialization (if already done)
;; allpoints -- ignore the transit epoch and consider all points -
;;              this is useful for KIC 1255 data
;; noLine - don't fit a line to the trend
;; psplot -- make a postscript plot

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/opt_correlations/opt_corr'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif


  if not keyword_set(fast) then begin
     get_profile_widths
     plot_tim_ser
  endif

  restore,'data/prof_widths.sav' ;; spatial widths
  restore,'data/timedata.sav' ;; time & orbital phase
  restore,'data/compiled_model_params.sav' ;; get the transit start & end

  ;; get the binned data
  restore,'data/specdata.sav'

  ;; Find the out of transit points
  if keyword_set(allpoints) then begin
     offp = lindgen(n_elements(tplot))
  endif else begin
     offp = where(tplot LT hstart OR tplot GT hend)
  endelse

  ;; Make a median spectrum to divide out & normalize time series
  nwavs = n_elements(bingrid)
  ntime = n_elements(tplot)
  meddivspec = fltarr(nwavs)
  for i=0l,nwavs-1l do begin
     meddivspec[i] = median(binfl[i,*])
  endfor
  replicatedspec = rebin(meddivspec,nwavs,ntime)
  binfl = binfl / replicatedspec


  for j=0,nwavbins-1l do begin
     wavname = string(bingrid[j],format='(F4.2)')+'-'+$
               string(bingrid[j]+binsizes[j],format='(F4.2)')+'um'
     ;; Ignore all points that are exactly 2.194
     goodp = where((widths[offp,1] LT 2.19 OR widths[offp,1] GE 2.20) and $
                  finite(binfl[j,offp]) and finite(binfl[j,offp]))
     xp = binfl[j,offp[goodp]]
     yp = widths[offp[goodp],1]

     plot,xp,yp,$
          xrange=threshold(xp,low=0.2,high=0.8),$
          yrange=threshold(yp,low=0.2,high=0.8),$
          xtitle='Flux Ratio',$
          ytitle='Spatial Width',$
          psym=4,title=wavname
     linefit = robust_linefit(xp,yp,yfit)
     if not keyword_set(noLine) then begin
        oplot,xp,yfit,color=mycol('yellow')
     endif

     wait,0.5

  endfor

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif


end
