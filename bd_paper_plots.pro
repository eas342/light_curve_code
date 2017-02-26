pro bd_paper_plots,type,psplot=psplot,manualShift=manualShift
;; Save plots for the BD paper
;; examples: specphot-0835
;; psplot - postscript output vs X11
;; manualShift - do manual wavelength shifting

if not keyword_set(psplot) then psplot=1
if not keyword_set(manualShift) then manualShift=0

if manualShift EQ 0 then specShift=1 else specShift=0

;; Change to the appropriate source and speclist
case 1 of
   strmatch(type,'*0835'): begin
      change_planets,pname='2massj0835'
      choose_speclist,fchoice='file_lists/2mass_j0835_es_red_01.txt'
      srcName = '2MASS J0835-0819'
   end
   strmatch(type,'*1821'): begin
      change_planets,pname='2massj1821'
      choose_speclist,fchoice='file_lists/2mass_1821_es_red01.txt'
      srcName = '2MASS J1821+1414'
   end
   else: begin
      print,"Unkown system"
      return
   end

endcase

case 1 of
   strmatch(type,'specphot*'): begin
      compile_spec,/readc,specsh=specShift,manualshift=manualShift
      double_specphot,targetStarName=srcName,/jd,/hr,psplot=psplot,$
                      specLabelCharsize=0.45,custfrange=[0,1.3]
   end
   strmatch(type,'tser*'): begin
      compile_spec,/readc,nwavbins=25,custrange=[0.87,2.38],$
                   specsh=specShift,manualshift=manualShift
      plot_tim_ser,/singlep,/jd,/lind,psplot=psplot,$
                   /tallplot,/littleCirc,/hr
   end
   strmatch(type,'saveclean*'): begin
      ;; Saves cleaned time series. Doesn't apply linear de-trending
      ;; I then pass this data to tser_tools for emcee fitting
      compile_spec,/readc,nwavbins=25,custrange=[0.87,2.38],$
                   manualshift=1
      plot_tim_ser,/singlep,/jd,$
                   /tallplot,/littleCirc,/hr
   end
   strmatch(type,'makeSpec*'): begin
      compile_spec,/readc,nwavbins=25,custrange=[0.87,2.38],$
                   specsh=specShift,manualshift=manualShift
      plot_tim_ser,/singlep,/jd,/lind,/tallplot,/littleCirc,/fitcurve,$
         /offtranserr,/sinfit
      plot_rad_vs_wavl,/amp,custxrange=[0.8,2.5],custyrange=[-0.1,1.5],$
         /shade,psplot=psplot
   end
   strmatch(type,'airmassCheck*'): begin
      compile_spec,/readc,nwavbins=1,custrange=[1.35,1.41],$
                   specsh=specShift,manualshift=manualShift
      plot_tim_ser,/jd,/singlep,/hr,/littleCirc,/fitcurve,/offtranserr,/sinfit,psplot=psplot
   end
   else: begin
      print,'Unknown BD paper plot'
      return
   end
endcase
end
