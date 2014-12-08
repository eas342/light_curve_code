pro multi_night_plots,psplot=psplot,$
                        custwavrange=custwavrange,$
                      starRatios=starRatios
;; Goes through all nights in order to plot the different things (such
;; as specphot images)
;; psplot - save a postcript plot for each night
;; custwavrange - pass on to compile_spec the custom wavelength range
;; starRatios - plot the ratios of the two stars

  case 1 of
     keyword_set(starRatios): begin
        plotprenm = 'plots/individual_spectra/all_indspec'
     end
     else: plotprenm = 'plots/specphot_images/all_specphots'
  endcase
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
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

   case 1 of
      keyword_set(starRatios): begin
         plot_stars,/divide,custyrange=[0.8,2.8],/nolegend,custtitle=usedate
      end
      else: begin
         plot_specphot,/useclean,/removel,/usebin,custtitle=usedate,$
                 custxmargin=[9,0]
      end
   endcase
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
     !p.thick=1
     !x.thick=1
     !y.thick=1
  endif
!p.multi=0
!x.omargin = [0,0]

end
