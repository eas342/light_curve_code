pro plot_rad_vs_wavl,psplot=psplot

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/rad_vs_wavl'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif


  ;; read in the radius versus wavelength file
  readcol,'radius_vs_wavelength/radius_vs_wavl.txt',wavl,rad,rade,skipline=1,format='(F,F,F)'

  plot,wavl,rad,$
       xtitle='Wavelength (um)',$
       ytitle='Rp/R*',$
;       ystyle=16,xstyle=1,$
       ystyle=16,xstyle=0,$
       yrange=[0.1,0.2]
  oploterror,wavl,rad,rade
  
  oploterror,wavl,rad,fltarr(n_elements(wavl)),rade ;errcolor=mycol('yellow'),$
;                color=mycol('yellow') 
  
  ;; As in Gibson et al. 2012, show 3 scale heights around the
  ;; adopted Rp/R* from Jacob Bean et al. 2012
  scaleH = 0.0023E
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433],color=mycol(['red'])
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]+3E*scaleH,color=mycol(['red']),linestyle=2
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]-3E*scaleH,color=mycol(['red']),linestyle=2
        
     
     

     if keyword_set(psplot) then begin
        device, /close
        cgPS2PDF,plotprenm+'.eps'
        spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
        device,decomposed=0
        set_plot,'x'
        !p.font=-1
     endif
end  
