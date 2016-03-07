pro kic_1255_paper_plots,type
;; This saves the routines for making plots for the KIC 1255 paper
;; type = 'dspecphot' double spectro-photometry
;; type = 'difftime' differential time series
case type of
   'aug18check': begin
      choose_speclist,fchoice='file_lists/speclist_kic1255_2014aug18_es_xtract.txt'
      compile_both,/readc,custrange=[1.48,1.78],nwavbins=1,/bothb
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],custsep=0.01,$
                   /fitcurve,/kepfit,/offtranserr,/lind,custyrange=[0.95,1.01]
      spawn,'open radius_vs_wavelength/radius_vs_wavl.txt'
   end
   'control': begin
      choose_speclist,fchoice='file_lists/speclist_kic1255_2014sep03es_xtract_short.txt'
      compile_both,/readc,nwavbins=5,/pretend,inject=0
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.02,/fitcurve,/kepfit,/offtranserr,/lind,$
                   deletePS=0
      kic1255_mult
   end
   'diffcontrol': begin
      choose_speclist,fchoice='file_lists/speclist_kic1255_2014sep03es_xtract_short.txt'
;      compile_both,/readc,nwavbins=5,/pretend,inject=0
      compile_spec,/readc,nwavbins=5,/pretend,inject=0
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.02,/fitcurve,/kepdiff,/offtranserr,/lind,$
                   deletePS=0,/diff
      kic1255_mult
   end
   'controlboot': begin
      choose_speclist,fchoice='file_lists/speclist_kic1255_2014sep03es_xtract_short.txt'
      compile_both,/readc,nwavbins=5,/pretend,inject=0
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.02,/fitcurve,/kepfit,/offtranserr,/lind,/boot
      spawn,'cp radius_vs_wavelength/radius_vs_wavl.txt radius_vs_wavelength/radius_vs_wavl_boot.txt'
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.02,/fitcurve,/kepfit,/offtranserr,/lind
      plot_rad_vs_wavl,/depthk,custyrange=[-0.5,1.4],totsets=2,custxmargin=[9,2],$
                       preset='param_input/rad_choices_bootst_vs_mpfit_covar.txt',/psplot,$
                       custxrange=[0.45,2.4]
   end
   'contdiff': begin
      ;; Finds a differential time series
      choose_speclist,fchoice='file_lists/speclist_kic1255_2014sep03es_xtract_short.txt'
      compile_spec,/readc,/specsh,nwavbins=5,/pretend ;; simulate as if fitting a transit
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.075,0.09],$
                   custsep=0.02,/fitcurve,/kepdiff,/offtranserr,/lind,/diff,$

      compile_both,/readc,nwavbins=5,/pretend,inject=0
;      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.08,0.095],$
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.02,/fitcurve,/kepfit,/offtranserr,/lind,$
                   deletePS=0
      kic1255_mult
   end
   'inddiff': begin
      plot_rad_vs_wavl,/depthk,custyrange=[-1.4,1.4],totsets=5,custxmargin=[9,2],/diff,$
                       preset='param_input/rad_choices_moris_big_diff.txt',/psplot,$
                       custxrange=[0.45,2.4]
;      plot_rad_vs_wavl,/depthk,custyrange=[-1,1],totsets=6,custxmargin=[9,2],/diff,$
;                       preset='param_input/rad_choices_6_diff.txt',/psplot,$
;                       custxrange=[0.45,2.4]
      spawn,'open -R plots/rad_vs_wavl.eps'
   end
   'inddiffboot': begin
      plot_rad_vs_wavl,/depthk,custyrange=[-1.4,1.4],totsets=5,custxmargin=[9,2],/diff,$
                       preset='param_input/rad_choices_moris_big_diff_boot.txt',/psplot,$
                       custxrange=[0.45,2.4]
      spawn,'open -R plots/rad_vs_wavl.eps'
   end
   'avgdiff': begin
      ;; Get the normalization from the average photometry
;      avg_radii,preset='param_input/rad_choices_6_absolute.txt',totsets=6
      avg_radii,preset='param_input/rad_choices_moris_big.txt',totsets=5
      readcol,'radius_vs_wavelength/avg_rp_rs.txt',skipline=1,format='(F,F,F,F)',$
              binmid,binwidth,rad,raderr
      norm = rad[0]
      readcol,'transit_info/kic1255_ksc_depth.txt',mult,format='(F)',skipline=1,/silent
      print,'MORIS constant= ',rad[0] * mult[0],' +/- ',raderr[0] * mult[0]
;      avg_radii,preset='param_input/rad_choices_6_diff.txt',totsets=6,diffconst=norm
      spawn,'cp radius_vs_wavelength/radius_vs_wavl.txt radius_vs_wavelength/radius_vs_wavl_boot.txt'
      avg_radii,preset='param_input/rad_choices_moris_big_diff.txt',totsets=5,diffconst=norm
      plot_rad_vs_wavl,/depthk,custyrange=[-0.2,0.8],custxmargin=[9,2],kepthick=4,$
                       preset='param_input/rad_choices_avg.txt',/showmie,/psplot
   end
   'avgdiffBoot': begin
      avg_radii,preset='param_input/rad_choices_moris_big.txt',totsets=5
      readcol,'radius_vs_wavelength/avg_rp_rs.txt',skipline=1,format='(F,F,F,F)',$
              binmid,binwidth,rad,raderr
      norm = rad[0]
      readcol,'transit_info/kic1255_ksc_depth.txt',mult,format='(F)',skipline=1,/silent
      print,'MORIS constant= ',rad[0] * mult[0],' +/- ',raderr[0] * mult[0]
      avg_radii,preset='param_input/rad_choices_moris_big_diff_boot.txt',totsets=5,diffconst=norm
      plot_rad_vs_wavl,/depthk,custyrange=[-0.2,0.8],custxmargin=[9,2],kepthick=4,$
                       preset='param_input/rad_choices_avg.txt',/showmie,/psplot
   end
   'avgdiffBootCompar': begin
      avg_radii,preset='param_input/rad_choices_moris_big.txt',totsets=5
      readcol,'radius_vs_wavelength/avg_rp_rs.txt',skipline=1,format='(F,F,F,F)',$
              binmid,binwidth,rad,raderr
      norm = rad[0]
      readcol,'transit_info/kic1255_ksc_depth.txt',mult,format='(F)',skipline=1,/silent
      print,'MORIS constant= ',rad[0] * mult[0],' +/- ',raderr[0] * mult[0]
      avg_radii,preset='param_input/rad_choices_moris_big_diff_boot.txt',totsets=5,diffconst=norm
      spawn,'cp radius_vs_wavelength/avg_rp_rs.txt radius_vs_wavelength/radius_vs_wavl_boot.txt'
      avg_radii,preset='param_input/rad_choices_moris_big_diff.txt',totsets=5,diffconst=norm
      spawn,'cp radius_vs_wavelength/avg_rp_rs.txt radius_vs_wavelength/radius_vs_wavl.txt'

      plot_rad_vs_wavl,/depthk,custyrange=[-0.2,0.7],totsets=2,custxmargin=[9,2],/diff,$
                       preset='param_input/rad_choices_bootst_vs_mpfit_covar.txt',/psplot,$
                       custxrange=[0.45,2.4]

;      plot_rad_vs_wavl,/depthk,custyrange=[-0.2,0.8],custxmargin=[9,2],kepthick=4,$
;                       preset='param_input/rad_choices_avg.txt',/showmie,/psplot
   end
   'photonly': begin
      spawn,'cp file_lists/multi_night_6night.txt file_lists/multi_night.txt'
      compile_multi_night,/photonly
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.12,0.13],custsep=0.01,/fitcurve,$
                   /kepfit,/offtranserr,/lind,custyrange=[0.996,1.003]
   end
   'photbig': begin
      spawn,'cp file_lists/multi_night_high_moris_nights.txt file_lists/multi_night.txt'
      compile_multi_night,/photonly
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.12,0.13],custsep=0.01,/fitcurve,$
                   /kepfit,/offtranserr,/lind,custyrange=[0.996,1.003],legord=0,/psplot,deletePS=0
   end
   'photall': begin
      spawn,'cp file_lists/multi_night_8night.txt file_lists/multi_night.txt'
      compile_multi_night,/photonly
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.12,0.13],custsep=0.01,/fitcurve,$
                   /kepfit,/offtranserr,/lind,custyrange=[0.993,1.003],legord=0,/psplot,deletePS=0,$
                   /showkep
   end
   '9moris': begin
      spawn,'cp file_lists/multi_night_9night.txt file_lists/multi_night.txt'
      multi_night_plots,/phot,/fixR,/psplot,/showbjd
   end
   '9spex': begin
      multi_night_plots,/stser,/psplot,/fixr,custwavrange=[0.82,2.4]
   end
   'difftime': begin
      spawn,'cp file_lists/multi_night_6night.txt file_lists/multi_night.txt'
      compile_multi_night,mnwavbins=5,/diff
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.01,/fitcurve,/kepdiff,/offtranserr,/lind,/diff,$
                   deletePS=0,/psplot
   end
   'bigmoris': begin
      ;spawn,'cp file_lists/multi_night_6night.txt file_lists/multi_night.txt'
      spawn,'cp file_lists/multi_night_high_moris_nights.txt file_lists/multi_night.txt'
      compile_multi_night,mnwavbins=5,/diff
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.01,/fitcurve,/kepdiff,/offtranserr,/lind,/diff,$
                   deletePS=0,/psplot
   end
   'absSpec': begin
      spawn,'cp file_lists/multi_night_6night.txt file_lists/multi_night.txt'
      compile_multi_night,mnwavbins=5
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.01,/fitcurve,/kepdiff,/offtranserr,/lind
      avg_radii,preset='param_input/rad_choices_6_absolute.txt',totsets=6
      plot_rad_vs_wavl,/depthk,custyrange=[-0.2,0.8],custxmargin=[9,2],kepthick=4,totsets=2,$
                       preset='param_input/rad_choices_two_avgs.txt',/psplot

   end
   'doubleSph': begin
      spawn,'cp file_lists/multi_night_6night.txt file_lists/multi_night.txt'
      set_specphot_range,0.98,1.005
      compile_multi_night,mnwavbins=25,/diff ;(only load spectroscopy)
      plot_tim_ser,/singlep,timebin=40,custxrange=[-0.11,0.11],/diff
      ;;(to choose a night where there is spectroscopy to generate the star
      ;;plots, I chose Sep 02, 2014)
      choose_speclist,fchoice='file_lists/speclist_kic1255_2014sep02es_xtract.txt'
      compile_spec,/readc,/specsh,nwavbins=25

      double_specphot,/useclean,/skipI,/psplot
      set_specphot_range,0.97,1.005
   end
   'oldcontrol': begin
      choose_speclist,fchoice='file_lists/speclist_kic1255_2014sep03es_xtract.txt'
      compile_both,nwavbins=5,/readC,/pretendTransit,inject=0
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],custsep=0.025,$
                   /fitcurve,/kepdiff,/offtranserr,/lind,/diff
      plot_rad_vs_wavl,/depthk,custyrange=[-0.5,0.5]
      kic1255_mult
   end
   'rms_bin_size': begin
      spawn,'cp file_lists/multi_night_9night.txt file_lists/multi_night.txt'
      multi_night_plots,/binsizephot,/fixr,/psplot
   end
   else: message,'Input not found.',/cont
endcase

end
