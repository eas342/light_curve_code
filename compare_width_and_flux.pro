pro compare_width_and_flux,fast=fast
;; Compares with spatial profile widths to flux
;; fast -- skips the initialization (if already done)

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
  offp = where(tplot LT hstart OR tplot GT hend)

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
     plot,binfl[j,offp],widths[offp,1],$
          xrange=[0.95,1.05],$
          yrange=[2,6],$
          xtitle='Flux Ratio',$
          ytitle='Spatial Width',$
          psym=4
     wait,0.5
  endfor

end
