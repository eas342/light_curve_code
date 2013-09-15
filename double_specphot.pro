pro double_specphot,psplot=psplot
;; Puts a plot of the stars spectra directly on top of the specphot
;; plot
;; psplot - saves a postscript plot


  plot_tim_ser

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotnmpre = 'plots/specphot_images/specphot_image'
     device,encapsulated=1, /helvetica,$
            filename=plotnmpre+'.eps',bits_per_pixel=8
     device,xsize=10, ysize=13,decomposed=1,/color
  endif

  !p.multi = [0,1,2]
  !p.position = [0.2,0.75,0.75,0.95]
  plot_stars,/normall,/showback,/directText,custXmargin=[9,12],/skipXTitle,$
             custYmargin=[0,0]
  !p.position = [0.2,0.1,0.75,0.748]
  plot_specphot,/removelin,/skipInitialize,custymargin=[4,4]
  !p.position = [0,0,0,0]
  !p.multi = 0

  if keyword_set(psplot) then begin
     device,/close
     cgPS2PDF,plotnmpre+'.eps'
     spawn,'convert -density 450% '+plotnmpre+'.pdf '+plotnmpre+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif


end
