pro kic_1255_paper_plots,type
;; This saves the routines for making plots for the KIC 1255 paper
;; type = 'dspecphot' double spectro-photometry
;; type = 'difftime' differential time series
case type of
   'difftime': begin
      spawn,'cp file_lists/multi_night_6night.txt file_lists/multi_night.txt'
      compile_multi_night,mnwavbins=5,/diff
      plot_tim_ser,timebin=40,/singlep,custxrange=[-0.2,0.13],$
                   custsep=0.01,/fitcurve,/kepdiff,/offtranserr,/lind,/diff,$
                   deletePS=0,/psplot
   end
   else: begin
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
endcase

end