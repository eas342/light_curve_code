pro plot_row_experiment,reinitialize=reinitialize,$
                        custyrange=custyrange,psplot=psplot
;;plots the row experiment
;; reinitialize - re-initializes the time series-gathering


  if keyword_set(psplot) then begin
     plotnm = 'plots/row_exp/row_exp_plot.eps'
     set_plot,'ps'
     !p.font=0
     !p.thick=3
     !x.thick=3
     !y.thick=3
     device,encapsulated=1, /helvetica,$
            filename=plotnm
     device,xsize=17,ysize=12,decomposed=1,/color
  endif

if keyword_set(reInitialize) then begin
   
   nwavbins=9
   ModArray = 0.04 * findgen(10) + 0.8
   nMods = n_elements(modArray)
   rmsArrExp = fltarr(nwavbins,nMods)
   
   for i=0l,nMods-1l do begin
      compile_spec,/readc,/specsh,specKey=7+i,nwavbins=nwavbins
      plot_tim_ser
      restore,'data/rmsdata.sav'
      rmsArrExp[*,i] = fracRMSarr
   endfor
   
   save,bingrid,fracRMSarr,tsizes,bingridmiddle,binsizes,$
        planetdat,fracPhotonarr,rmsArrExp,modArray,$
        filename='data/rmsdata_rowexp.sav'
endif else begin
   restore,'data/rmsdata_rowexp.sav'
endelse

nMods = n_elements(modArray)

if n_elements(custyrange) EQ 0 then custyrange=[min(rmsArrexp),max(rmsArrexp)]*100E

plot,modArray,rmsArrExp[0,*]*100E,$
     yrange=custyrange,$
     xtitle='Row Mult factor',ytitle='Robust RMS (%)',$
     /nodata
arraycol = myarraycol(nMods,psversion=psplot)
nwavs = n_elements(bingrid)
for i=0l,nwavs-1l do begin
   oplot,modArray,rmsArrExp[i,*]*100E,color=arraycol[i]
   polyTrend = poly_fit(modArray,transpose(rmsArrExp[i,*]),2)
   trough = -polytrend[1]/(2E * polyTrend[2])
;   minRMS = min(rmsArrExp[i,*],minind)
   oplot,[trough],[poly(trough,polyTrend)]*100E,psym=4
   
endfor
al_legend,string(bingrid,format='(F6.2)')+'um',color=arraycol,linestyle=0,$
          /clear

  if keyword_set(psplot) then begin
     device,/close
     ;; return the plotting to normal
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick=1
     !x.thick=1
     !y.thick=1
  endif

end
