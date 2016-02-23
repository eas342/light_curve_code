pro multi_night_plots,psplot=psplot,$
                        custwavrange=custwavrange,$
                      starRatios=starRatios,fitdepths=fitdepths,$
                      statep=statep,fixrange=fixrange,$
                      indWav=indWav,photometry=photometry,$
                      differential=differential,$
                      starplots=starplots,indflux=indflux,$
                      nightsummary=nightsummary,boot=boot,stser=stser,$
                      binsizephot=binsizephot
;; Goes through all nights in order to plot the different things (such
;; as specphot images)
;; psplot - save a postcript plot for each night
;; custwavrange - pass on to compile_spec the custom wavelength range
;; starRatios - plot the ratios of the two stars
;; fitdepths - fit the transit depths and save all
;; statep - save the state parameter plots
;; indWav - show the individual wavelengths
;; photometry - just show the photometry
;; differential - find the differential spectrum
;; indflux - show the individual flux
;; nightsummary - print out a summary of the nights
;; boot - use bootstrap errors
;; stser - SpeX time series
;; binsizephot - show the RMS versus bin size

  case 1 of
     keyword_set(stser): begin
        plotprenm = 'plots/spec_t_series/all_spex'
     end
     keyword_set(photometry): begin
        plotprenm = 'plots/spec_t_series/all_phot'
     end
     keyword_set(indflux): begin
        plotprenm = 'plots/spec_t_series/individ_flux'
     end
     keyword_set(starplots): begin
        plotprenm = 'plots/individual_spectra/all_indspec'
     end
     keyword_set(starRatios): begin
        plotprenm = 'plots/individual_spectra/all_indspec_ratio'
     end
     keyword_set(binsizephot): begin
        plotprenm = 'plots/binsize/all_binsize'
     end
     else: plotprenm = 'plots/specphot_images/all_specphots'
  endcase
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps',bits_per_pixel=8
;      device,xsize=18, ysize=10,decomposed=1,/color
      device,xsize=18, ysize=20,decomposed=1,/color
     !p.thick=3
     !x.thick=3
     !y.thick=3
  endif

if keyword_set(indWav) then usebin=0 else usebin=1

;; get the list of speclists
  readcol,'file_lists/multi_night.txt',listFile,format='(A)',comment='#'

nNights = n_elements(listFile)

!p.multi = [0,3,3]
if keyword_set(photometry) then begin
   !x.omargin=[0,0]
endif else !x.omargin = [0,20]

for i=0l,nNights-1l do begin
   choose_speclist,fchoice=listFile[i]

   ;; Get the spectral data
   if n_elements(custwavrange) EQ 0 then begin
      if keyword_set(statep) then begin
         custwavrange = [1.1,1.4]
      endif else custwavrange = [0.95,2.35]
   endif

   restore,'data/used_date.sav'

   case 1 of
      keyword_set(photometry): begin
         if strmatch(usedate,'*2014aug18*') then begin
            compile_phot,/readc,/both
         endif else compile_phot,/readC
      end
      keyword_set(stser): begin
         compile_spec,/readC,nwavbins=5,custrange=custwavrange,/specsh
      end
      keyword_set(statep): begin
         get_profile_widths,/esX
         compile_spec,/readC,nwavbins=1,custrange=custwavrange,/specshift,/saveshifts,/quickread
      end
      keyword_set(indflux): begin
         compile_spec,nwavbins=5,/readc,/quickread,/specsh;,custrange=[0.9,1.2]
      end
      keyword_set(binsizephot): begin
         compile_phot,/readc
      end
      else: begin
         compile_spec,/specsh,/readC,custrange=custwavrange,nwavbins=25,/quickread
      end
   endcase

   if strmatch(usedate,'*2014sep03*') then secondary=1 else secondary=0
   if keyword_set(fixRange) then begin
      case 1 of
         secondary EQ 1: begin
            custSpecPhRange = [.37,.55]
         end
         else: custSpecPhRange = [-0.1,0.13]
      endcase
   endif else undefine,custSpecPhRange

   plot_tim_ser,timebin=40,/lind,/offtranserr,/noplots,secondary=secondary,$
                custXrange=custSpecPhRange

   case 1 of
      keyword_set(photometry): begin
         if showdate EQ '2014 Sep03' then showdate = showdate +' (Control)'
         if showdate EQ '2013 Aug13' then labelKep=1 else labelKep=0
         plot_tim_ser,timebin=40,/lind,/offtranserr,secondary=secondary,$
                      custXrange=custSpecPhRange,/showkep,custtitle=showdate,$
                      custyrange=[0.983,1.0045],/singlep,/skipreset,custsep=0.0075,$
                      labelKep=labelKep
         ;; use the same variable as custom specphot range
         
      end
      keyword_set(stser): begin
         if i EQ 3 - 1 then skipwavl=0 else skipwavl=1
         plot_tim_ser,timebin=40,/lind,/offtranserr,secondary=secondary,skipwavl=skipwavl,$
                      custXrange=custSpecPhRange,/singlep,/skipreset,custtitle=showdate
      end
      keyword_set(indflux): begin
         plot_tim_ser,/ind,custtitle=showdate,secondary=secondary,$
                      /singlep,/skipreset,custsep=0.05 ;custyrange=[0.9,1.15],
      end
      keyword_set(starRatios): begin
         plot_stars,/divide,custyrange=[0.8,2.8],/nolegend,custtitle=showdate
      end
      keyword_set(starplots): begin
         plot_stars,/nolegend,custtitle=showdate,/showb
      end
      keyword_set(statep): begin
         state_parameters,/psplot
         spawn,'cp plots/spec_t_series/state_param.eps '+$
               'plots/state_params/state_p_'+usedate+'.eps'
      end
      keyword_set(fitDepths): begin
         compile_both,nwavbins=5,/readc,/specsh
         plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                      custsep=0.01,/fitcurve,/kepfit,/offtranserr,/lind,legord=1
         spawn,'cp radius_vs_wavelength/radius_vs_wavl.txt '+$
               'radius_vs_wavelength/rad_vs_wavl_'+usedate+'.txt'
      end
      keyword_set(differential): begin
         compile_spec,nwavbins=5,/readc,/specsh,/quickread
         plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                      custsep=0.01,/fitcurve,/kepdiff,/offtranserr,/lind,/diff,legord=1,boot=boot
         if keyword_set(boot) then addDescrip='boot_' else addDescrip=''
         spawn,'cp radius_vs_wavelength/radius_vs_wavl.txt '+$
               'radius_vs_wavelength/'+addDescrip+'diff_spec_'+usedate+'.txt'
      end
      keyword_set(nightsummary): begin
         get_profile_widths,/useSaved
         restore,'data/prof_widths.sav'
         restore,'data/specdata.sav'
         if median(utgrid) LT date_conv('2014-07-25T12:00:00.00','JD') then begin
            ;; Pre-upgrade for detector
            PSDetect = 0.15E
         endif else begin
            PSDetect = 0.1E
         endelse
         medSeeing = median(widths * 2.35E * PSdetect) ;; FWHM in arcsec
         medExpTime = median(itimeGrid) ;; in sec
         restore,'data/timedata.sav'
         tingress = (hstart - min(tplot)) * planetdat.period * 24D
         tegress = (max(tplot) - hend) * planetdat.period * 24D

         if i EQ 0 then begin
            openw,1,'data/night_summary.txt'
            printf,1,'Date','FWHM','T_exp','T_ingress','T_egress',$
                   'Airmass Range','Sun Alt',$
                   format='(A12,4A10,A15,A10)'
         endif
         lastair = airmass[n_elements(airmass)-1l]
         airstring = string(airmass[0],min(airmass),lastair,format='(F4.2,"-",F4.2,"-",F4.2)')
         sunstring = string(altsun[0],altsun[n_elements(altsun)-1l],$
                            format='(I3," - ",I3)')
         printf,1,useDate,medSeeing,medExpTime,tingress,tegress,$
                airstring,sunstring,$
                format='(A12,F10.2,3F10.1,A15,A10)'
      end
      keyword_set(binsizephot): begin
         plot_bin_size,/photmode,custyrange=[1E-2,1],custtitle=showdate,$
                       tSerRange=custSpecPhRange,/nointerp
      end
      else: begin
         plot_specphot,usebin=usebin,/removel,custtitle=showdate,$
                       custxmargin=[9,0],/skipI,secondary=secondary,$
                       custyrange=custSpecPhRange,thickmarkers=2
      end
   endcase

;,/psplot

;   spawn,'mv plots/specphot_images/specphot_image.png plots/specphot_images/specphot_image_'+usedate+'.png'

endfor
close,1
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
