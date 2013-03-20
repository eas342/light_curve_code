pro plot_multiple_epochs

  ;; get the compiled spectroscopic data
  restore,'data/specdata.sav'

  ;; get the time info
  restore,'data/timedata.sav'

  IndArr = [0,100,200,300,400,500]
;  IndArr = [0,300]
;  colorarr=mycol(['white','yellow'])
  colorarr=mycol(['white','yellow','green','lblue','orange','red'])

  tempXrange=[0.8,2.5]
  plot,lamgrid,flgrid[*,0,IndArr[0]],ytitle='Flux DN',charsize=2,$
       font=1,xtitle='Wavelength (um)',$
       /nodata,xmargin=[8,2],xstyle=1,xrange=tempXrange,$
       yrange=[0,1.5e6]

  for i=0l,n_elements(IndArr)-1l do begin
     oplot,lamgrid,flgrid[*,1,Indarr[i]],color=colorarr[i],psym=10
     oplot,lamgrid,flgrid[*,0,Indarr[i]],color=colorarr[i],psym=10
  endfor

  legend,'Orbital Phase = '+string(tplot[IndArr]),color=colorarr,$
         linestyle=0+intarr(n_elements(IndArr)),/right

end
