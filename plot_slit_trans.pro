pro plot_slit_trans,psplot=psplot,voigt=voigt,closerparams=closerparams
;; Shows the slit transmission as a function of other parameters
;; psplot - saves postscript & png plots
;; voigt - use a voigt profile with a damping parameter a=voigt instead of Gaussian
;; closerparams - parameters close to the measured ones

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/slit_trans/general_trans'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     device,xsize=15, ysize=10,decomposed=1,/color
     !p.thick=2
     !x.thick=2
     !y.thick=2
  endif

  ;; Make array of positions
  npts = 64
  H = 10E ;; half slit width
  if keyword_set(closerParams) then begin
     sigma = [1E,2E,3E,4E]
     y = dindgen(npts)/double(npts) * 6E - 3E
     myYrange = [80,97]
  endif else begin
     y = dindgen(npts)/double(npts) * 30E - 15E
     sigma = [3E,5E,10E,15E]
     myYrange=[0,0]
  endelse
  sigmaString = cgGreek('sigma')+' = '+string(sigma,format='(F8.1)')+' px '
  nlines = n_elements(sigma)
  farray = dblarr(nlines,npts)
  for i=0l,nlines-1l do begin
     if keyword_set(voigt) then begin
        farray[i,*] = voigt_slit(y,H,sigma[i],voigt)
     endif else farray[i,*] = gauss_slit(y,H,sigma[i])
  endfor
  
  colorArray = [!p.color,mycol(['red','blue','dgreen'])]
  linestyleArray = [0,1,3,4]
  plot,y,farray[0,*] * 100E,ytitle='Transmission (%)',$
       xtitle='Position (pixels)',$
       ystyle=16,/nodata,yrange=myYrange
  for i=0l,nlines-1l do begin
     oplot,y,farray[i,*]*100E,color=colorArray[i],linestyle=lineStyleArray[i]
  endfor

  if keyword_set(psplot) then myCharsize=0.7 else myCharsize=1.0

  al_legend,sigmastring,color=colorArray,linestyle=linestyleArray,$
            charsize=myCharsize
  if keyword_set(voigt) then begin
     al_legend,['H = 10px','a = '+string(voigt[0],format='(F6.2)')],$
               /right,charsize=myCharsize
  endif else begin
     al_legend,['H = 10px'],/right,charsize=myCharsize
  endelse

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 300% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick=1
     !x.thick=1
     !y.thick=1
  endif


end

