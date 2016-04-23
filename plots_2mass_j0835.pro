pro plots_2mass_j0835,type

if n_elements(type) EQ 0 then type='time series sine fit 01'
case type of
   'broadband sinefit': begin
      ;; This is the plot I sent out in Quick Look notes 02
      choose_speclist,fchoice='file_lists/2mass_j0835_es_red_01.txt'
      compile_spec,/readc,nwavbins=1,custrange=[1.0,2.3],/specsh
      plot_tim_ser,/singlep,/lind,timebin=40,custsep=0.01,/offtranserr,/sinfit,/fitcurve,custyrange=[0.992,1.007],/psplot
   end
   else: begin
      ;; This is the plot I sent out in Quick Look notes 02
      choose_speclist,fchoice='file_lists/2mass_j0835_es_red_01.txt'
      compile_spec,/readc,nwavbins=5,custrange=[1.0,2.3],/specsh
      plot_tim_ser,/singlep,/lind,timebin=40,custsep=0.01,/offtranserr,/sinfit,/fitcurve,/psplot
   end
endcase

end
