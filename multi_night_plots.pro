pro multi_night_specphot,psplot=psplot,$
                        custwavrange=custwavrange
;; Goes through all nights in order to plot the correlation between
;; position and flux ratio
;; psplot - save a postcript plot for each night
;; custwavrange - pass on to compile_spec the custom wavelength range


  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/specphot_images/all_specphots'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
      device,xsize=18, ysize=10,decomposed=1,/color
     !p.thick=3
     !x.thick=3
     !y.thick=3
  endif

;; get the list of speclists
  readcol,'file_lists/multi_night.txt',listFile,format='(A)',comment='#'

nNights = n_elements(listFile)

!p.multi = [0,4,2]
!x.omargin = [0,20]

for i=0l,nNights-1l do begin
   choose_speclist,fchoice=listFile[i]

   ;; Get the spectral data
   if n_elements(custwavrange) EQ 0 then custwavrange = [0.95,2.35]
   compile_spec,/readC,/specsh,custrange=custwavrange,nwavbins=25
   plot_tim_ser,timebin=40,/lind,/offtranserr,/noplots

   restore,'data/used_date.sav'

   plot_specphot,/useclean,/removel,/usebin,custtitle=usedate,$
                 custxmargin=[9,0]
;,/psplot

;   spawn,'mv plots/specphot_images/specphot_image.png plots/specphot_images/specphot_image_'+usedate+'.png'

endfor

  if keyword_set(psplot) then begin
     device,/close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 450% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
!p.multi=0
!x.omargin = [0,0]

end
