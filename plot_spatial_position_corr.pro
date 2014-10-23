pro plot_spatial_position_corr,psplot=psplot,$
                        custwavrange=custwavrange
;; Goes through all nights in order to plot the correlation between
;; position and flux ratio
;; psplot - save a postcript plot for each night
;; custwavrange - pass on to compile_spec the custom wavelength range

;; get the list of speclists
  readcol,'file_lists/multi_night.txt',listFile,format='(A)',comment='#'

nNights = n_elements(listFile)

for i=0l,nNights-1l do begin
   choose_speclist,fchoice=listFile[i]

   ;; Get the spectral data
   if n_elements(custwavrange) EQ 0 then custwavrange = [0.95,1.4]
   compile_spec,/readC,/specsh,custrange=custwavrange,nwavbins=1
   plot_tim_ser
   state_parameters

   restore,'data/used_date.sav'
   restore,'data/state_parameters/full_parameters/'+specfileListNamePrefix+'.sav'

;   restore,'data/specdata.sav'

   x = statepstruct.position1
   y = statepstruct.fluxratio
   

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/spatial_pos_corr/'+usedate+'_sp_corr'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     device,xsize=15, ysize=10,decomposed=1,/color
     !p.thick=3
     !x.thick=3
     !y.thick=3
  endif

   plot,x,y,psym=4,$
        xrange=threshold(x,low=0.1,high=0.9),$
        yrange=threshold(y,low=0.1,high=0.9),$
        xtitle='Spatial Position',ytitle='Flux Ratio'
   al_legend,[usedate]
  if keyword_set(psplot) then begin
     device, /close
;     cgPS2PDF,plotprenm+'.eps'
;     spawn,'convert -density 300% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick=1
     !x.thick=1
     !y.thick=1
  endif


endfor
end
