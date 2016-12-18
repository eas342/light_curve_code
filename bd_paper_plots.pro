pro bd_paper_plots,type,psplot=psplot
;; Save plots for the BD paper

if not keyword_set(psplot) then psplot=1
case type of
   'specphot0835': begin
      change_planets,pname='2massj0835'
      choose_speclist,fchoice='file_lists/2mass_j0835_es_red_01.txt'
      compile_spec,/readc,/specsh
      double_specphot,targetStarName='2MASS J0835-0819',/jd,/hr,psplot=psplot,$
                      specLabelCharsize=0.45,custfrange=[0,1.3]
   end
   'specphot1821': begin
      change_planets,pname='2massj1821'
      choose_speclist,fchoice='file_lists/2mass_1821_es_red01.txt'
      compile_spec,/readc,/specsh
      double_specphot,targetStarName='2MASS J1821+1414',/jd,/hr,psplot=psplot,$
                      specLabelCharsize=0.45,custfrange=[0,1.3]
   end
   'tser0835': begin
      change_planets,pname='2massj0835'
      choose_speclist,fchoice='file_lists/2mass_j0835_es_red_01.txt'
      compile_spec,/readc,/specsh,nwavbins=25,custrange=[0.87,2.38]
      plot_tim_ser,/singlep,/jd,/lind,psplot=psplot,$
                   /tallplot,/littleCirc,/hr
   end
   'tser1821': begin
      change_planets,pname='2massj1821'
      choose_speclist,fchoice='file_lists/2mass_1821_es_red01.txt'
      compile_spec,/readc,/specsh,nwavbins=25,custrange=[0.87,2.38]
      plot_tim_ser,/singlep,/jd,/lind,psplot=psplot,$
                   /tallplot,/littleCirc,/hr
   end
   'makeSpec1821': begin
      change_planets,pname='2massj1821'
      choose_speclist,fchoice='file_lists/2mass_1821_es_red01.txt'
      compile_spec,/readc,/specsh,nwavbins=25,custrange=[0.87,2.38]
      plot_tim_ser,/singlep,/jd,/lind,/tallplot,/littleCirc,/fitcurve,$
         /offtranserr,/sinfit
      plot_rad_vs_wavl,/amp,custxrange=[0.8,2.5],custyrange=[-0.1,1.5],$
         /shade,psplot=psplot
   end
   else: begin
      print,'Unknown BD paper plot'
   end
endcase
end
