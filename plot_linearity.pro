pro plot_linearity,psplot=psplot
;; the psplot keyword tells the program to make a Postscript file
  restore,'data/fluxsave.sav'
  ;; should pick up flarr(average pixel flux per file and box),
  ;;t(time),expt(seconds),nbox(number of
  ;; extraction boxes) and nfile(number of files),flerr(error in
  ;; average pixel flux),wavl(the approximate wavelength for each box)
  
  ;; convert expt from a string to a float
  expt = float(expt)

  if keyword_set(psplot) then begin
     plotnm = 'plots/basic_linearity_plot.eps'
     set_plot,'ps'
     !p.font=0
     device,encapsulated=1, /helvetica,$
            filename=plotnm
     device,xsize=14,ysize=10,decomposed=1,/color
  endif

  ylabel = 'Counts (DN)'
  y = flarr

  cgplot,expt,y[*,0],psym=2,$
         xtitle='Exposure Time (s)',$
         ytitle='Counts (DN)',/nodata,$
         xrange=[0,110]

  ;; make arrays for symbols and colors
  psymarr = symbol_colors(nbox,type='psym')
  lstylearr = symbol_colors(nbox,type='linestyle')
  colorarr = symbol_colors(nbox,type='color')

  ;; plot points for the line
  ppts = 20
  maxt = max(expt)
  mint = min(expt)
  plotx = indgen(ppts)/float(ppts-1)*(maxt)

  linefitarr = fltarr(nbox,2)
  linefitarre = fltarr(nbox,2) ;; errors in line fit parameters
  for i=0l,nbox-1l do begin
     oplot,expt,y[*,i],psym=psymarr[i],color=colorarr[i]
     ;; fit to a line
     linefitarr[i,*] = poly_fit(float(expt),flarr[*,i],1)
     ;; plot the line
     ploty = linefitarr[i,0] + plotx*linefitarr[i,1]
     oplot,plotx,ploty,linestyle=lstylearr[i],color=colorarr[i]
  endfor
 ;; Make a box legend
  if keyword_set(psplot) then begin
     plotnm2 = 'plots/legend_plot.eps'
     device,/close
     device,filenam=plotnm2
  endif
  cgplot,indgen(10),/nodata
  legend,'Box '+strtrim(indgen(nbox),1),$
         linestyle=lstylearr,$
         color=colorarr
  legend,'Box '+strtrim(indgen(nbox),1),$
         psym=psymarr,$
         color=colorarr,/right

  if keyword_set(psplot) then begin
     plotnm3 = 'plots/normalized_response_plot.eps'
     device,/close
     device,filenam=plotnm3
  endif

  ;; Plot the DN normalized by t_100
  cgplot,expt/expt[nfile-1l],y[*,0]/y[nfile-1l,0],psym=2,$
         xtitle='t / t!D 100 !N',$
         ytitle='(DN - bias) / DN!D 100 !N',$
         /nodata,$
         xrange=[0,1.1],yrange=[0,1.1]
  
  plotx = indgen(ppts)/float(ppts-1)*(1.1)

  ;; find the points where the exposure time is 100sec
  normpt = where(expt EQ 100.0)

  ;; make arrays for the average flux in the 100 second bin for each box
  avgfl = fltarr(nbox)
  avgfle = fltarr(nbox)

  ;; make arrays for a linear fit to the response
  respfitarr = fltarr(nbox,2)
  respfitarre = fltarr(nbox,2) ;; errors in line fit parameters


  for i=0,nbox-1l do begin
     x = expt/expt[nfile-1l]
     avgfl[i] = mean(flarr[normpt,i])
     ;; add the errors in quadrature
     avgfle[i] = sqrt(total(flerr[normpt,i]^2)) / float(n_elements(normpt))
     ;; normalize the flux
     y = (flarr[*,i] - linefitarr[i,0]) / avgfl[i]
     ;; add the fractioanl errors in quadrature
     yerr = sqrt((flerr[*,i]/flarr[*,i])^2 + (avgfle[i]/avgfl[i])^2) * y

     oplot,x,y,psym=psymarr[i],$
           color=colorarr[i]
     ;; fit to a line
     respfitarr[i,*] = poly_fit(x,y,1,sigma=paramsig,$
                               measure_errors=yerr)
     ;; plot the line
     ploty = respfitarr[i,0] + plotx*respfitarr[i,1]
     oplot,plotx,ploty,linestyle=lstylearr[i],color=colorarr[i]

     respfitarre[i,*] = paramsig

  endfor

  if keyword_set(psplot) then begin
     plotnm4 = 'plots/response_vs_dn.eps'
     device,/close
     device,filename=plotnm4
  endif

  ;; plot the slope versus the average counts to see if there is a
  ;; dependence
  cgplot,avgfl,respfitarr[*,1],$
         xtitle='DN!D 100 !N',$
         ytitle='Response Slope',$
         ystyle=16,/nodata,$
         xrange=[0,4000]
  oploterror,avgfl,respfitarr[*,1],fltarr(nbox),respfitarre[*,1],psym=2

  ;; print out the STDEV in response slope
  print,'Stdev in response slope: ',stddev(respfitarr[*,1])

  ;; plot the response as a function of wavelength to find any patterns
  if keyword_set(psplot) then begin
     device,/close
     plotnm5 = 'plots/slope_response_vs_wavl.eps'
     device,filename=plotnm5
  endif
  cgplot,wavl,respfitarr[*,1],$
         xtitle='Wavl (um)',$
         ytitle='Response Slope',$
         ystyle=16,/nodata
  oploterror,wavl,respfitarr[*,1],fltarr(nbox),respfitarre[*,1],psym=2

  if keyword_set(psplot) then begin
     device,/close
     cgPS2PDF,plotnm
     cgPS2PDF,plotnm2
     cgPS2PDF,plotnm3
     cgPS2PDF,plotnm4
     cgPS2PDF,plotnm5
     ;; return the plotting to normal
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

end
