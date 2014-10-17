pro plot_row_experiment,reinitialize=reinitialize
;;plots the row experiment
;; reinitialize - re-initializes the time series-gathering

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
plot,modArray,rmsArrExp[0,*]*100E,$
     yrange=[min(rmsArrexp),max(rmsArrexp)]*100E,$
     xtitle='Row Mult factor',ytitle='Robust RMS (%)',$
     /nodata
arraycol = myarraycol(nMods)
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

end
