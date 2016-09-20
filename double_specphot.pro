pro double_specphot,psplot=psplot,noremovelin=noremovelin,$
                    custxrange=custxrange,useclean=useclean,$
                    showmod=showmod,skipInitialize=skipInitialize,$
                    uniformPxGrid=uniformPxgrid,targetStarName=targetStarName
;; Puts a plot of the stars spectra directly on top of the specphot
;; plot
;; psplot - saves a postscript plot
;; noremovelin - do not perform a linear de-trend of the baseline
;; custxrange -- pass the custom xrange to plot_specphot
;; useclean - passed onto plot_specphot
;; showmod - passed onto plot_specphot
;; uniformPxgrid - force a uniform pixel grid (or a distorted wavelength grid)

if keyword_set(noremovelin) then begin
   myRemovelin = 0
endif else myRemovelin=1

if n_elements(custxrange) EQ 0 then custXrange=[0.95,2.4]

  if not keyword_set(useclean) then plot_tim_ser,/noplots

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotnmpre = 'plots/specphot_images/specphot_image'
     device,encapsulated=1, /helvetica,$
            filename=plotnmpre+'.eps',bits_per_pixel=8
     device,xsize=10, ysize=13,decomposed=1,/color
     !p.thick=3
  endif

  !p.multi = [0,1,2]
  !p.position = [0.2,0.75,0.75,0.95]
  plot_stars,/normall,/showback,/directText,custXmargin=[9,12],/skipXTitle,$
             custYmargin=[0,0],custxrange=custxrange,uniformPxGrid=uniformPxGrid,$
             targetStarName=targetStarName
  !p.position = [0.2,0.1,0.75,0.748]
  plot_specphot,removelin=myRemovelin,/skipInitialize,custymargin=[4,4],$
                custxrange=custxrange,useclean=useclean,showmod=showmod,$
                uniformWgrid=1
  !p.position = [0,0,0,0]
  !p.multi = 0

  if keyword_set(psplot) then begin
     device,/close
     cgPS2PDF,plotnmpre+'.eps'
     spawn,'convert -density 450% '+plotnmpre+'.pdf '+plotnmpre+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick=1
  endif


end
