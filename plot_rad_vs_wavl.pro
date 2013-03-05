pro plot_rad_vs_wavl,psplot=psplot,showstarspec=showstarspec
;;psplot -- saves a postscript plot
;;showstarspec -- shows a star spectrum on the same plot

  !x.margin = [13,14]
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

  ;; restore the star spectrum to show where telluric & star features
  ;; might be
  restore,'data/specdata.sav'

  if keyword_set(showstarspec) then ytempstyle=8 else ytempstyle=1

  plot,wavl,rad,$
       xtitle='Wavelength (um)',$
       ytitle='Rp/R*',$
;       ystyle=16,xstyle=1,$
       ystyle=ytempstyle,xstyle=1,xrange=[0.8,2.55],$
       yrange=[0.10,0.20],/nodata
;       yrange=[0.12,0.16],/nodata
;  oploterror,wavl,rad,rade
  wavlwidth = binsizes/2E
  oploterror,wavl,rad,wavlwidth,rade,psym=3
;                color=mycol('yellow') 
  

  

  ;; As in Gibson et al. 2012, show 3 scale heights around the
  ;; adopted Rp/R* from Jacob Bean et al. 2012
  scaleH = 0.00115E
;  scaleH = 0.00115E * 2E
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433],color=mycol(['red'])
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]+3E*scaleH,color=mycol(['red']),linestyle=2
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]-3E*scaleH,color=mycol(['red']),linestyle=2
        
     
     
  prevXrange = !x.crange
  if keyword_set(showstarspec) then begin
     ;; plot the source spectrum
     plot,lamgrid,flgrid(*,0,1),/noerase,xrange=prevXrange,ystyle=5,xstyle=1,$
          yrange=[-6E5,6E5],/nodata
     oplot,lamgrid,flgrid(*,0,1),color=mycol('blue')
     axis,yaxis=1,yrange=!y.crange,color=mycol('blue'),/ystyle,$
          ytitle='Raw Source Flux (DN)'
     
     !x.margin = [10.0,3.0]
  endif
  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
  
end  
  
