pro multi_night_plots,psplot=psplot,$
                        custwavrange=custwavrange,$
                      starRatios=starRatios,fitdepths=fitdepths,$
                      statep=statep
;; Goes through all nights in order to plot the different things (such
;; as specphot images)
;; psplot - save a postcript plot for each night
;; custwavrange - pass on to compile_spec the custom wavelength range
;; starRatios - plot the ratios of the two stars
;; fitdepths - fit the transit depths and save all
;; statep - save the state parameter plots

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
;      device,xsize=18, ysize=10,decomposed=1,/color
      device,xsize=18, ysize=20,decomposed=1,/color
     !p.thick=3
     !x.thick=3
     !y.thick=3
  endif

;; get the list of speclists
  readcol,'file_lists/multi_night.txt',listFile,format='(A)',comment='#'

nNights = n_elements(listFile)

!p.multi = [0,3,3]
!x.omargin = [0,20]

for i=0l,nNights-1l do begin
   choose_speclist,fchoice=listFile[i]

   ;; Get the spectral data
   if n_elements(custwavrange) EQ 0 then begin
      if keyword_set(statep) then begin
         custwavrange = [1.1,1.4]
      endif else custwavrange = [0.95,2.35]
   endif

   restore,'data/used_date.sav'

   if keyword_set(statep) then begin
      get_profile_widths,/esX
      compile_spec,/readC,nwavbins=1,custrange=custwavrange,/specshift,/saveshifts
   endif else begin
      compile_spec,/readC,/specsh,custrange=custwavrange,nwavbins=25
   endelse
   if strmatch(usedate,'*2014sep03*') then secondary=1 else secondary=0
   plot_tim_ser,timebin=40,/lind,/offtranserr,/noplots,secondary=secondary



   case 1 of
      keyword_set(starRatios): begin
         plot_stars,/divide,custyrange=[0.8,2.8],/nolegend,custtitle=usedate
      end
      keyword_set(statep): begin
         state_parameters,/psplot
         spawn,'cp plots/spec_t_series/state_param.eps '+$
               'plots/state_params/state_p_'+usedate+'.eps'
      end
      keyword_set(fitDepths): begin
         compile_both,nwavbins=5,/readc
         plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                      custsep=0.01,/fitcurve,/kepfit,/offtranserr,/lind,legord=2
         spawn,'cp radius_vs_wavelength/radius_vs_wavl.txt '+$
               'radius_vs_wavelength/radius_vs_wavl.txt/rad_vs_wavl_'+usedate+'.txt'
      end
      else: begin
         plot_specphot,/usebin,/removel,custtitle=usedate,$
                 custxmargin=[9,0],/skipI,secondary=secondary
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
