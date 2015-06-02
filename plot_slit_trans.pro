pro plot_slit_trans,psplot=psplot,voigt=voigt,$
                    closerparams=closerparams,$
                    avoigt=avoigt,differential=differential,$
                    starsep=starsep,FWHM=FWHM
;; Shows the slit transmission as a function of other parameters
;; psplot - saves postscript & png plots
;; voigt - use a voigt profile with a damping parameter a=voigt instead of Gaussian
;; closerparams - parameters close to the measured ones
;; avoigt -- use the simple analytic approximation to the Voigt function
;; differential  - show a two star differential (b/c that's what matters)
;; starsep - set the number of pixels the two stars are offset
;; FWHM -- show the FWHM instead of sigma

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/slit_trans/general_trans'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     device,xsize=15, ysize=10,decomposed=1,/color
     !p.thick=3
     !x.thick=3
     !y.thick=3
  endif

  ;; Make array of positions
  npts = 64
  H = 10E ;; half slit width
  if n_elements(starSep) EQ 0 then starSep = 1E ;; default separation between the two stars

  case 1 of
     keyword_set(voigt): vdamp = voigt
     keyword_set(avoigt): vdamp = avoigt
     else: vdamp = 0.0E
  endcase

  case 1 of
     keyword_set(closerParams): begin
;        sigma = [0.7E,1E,1.5E,2E]
        sigma = [3E,4E,5E,7E]/2.35E
        y = dindgen(npts)/double(npts) * 6E - 3E
        myYrange = [73,96]
     end
     else: begin
        y = dindgen(npts)/double(npts) * 30E - 15E
        sigma = [3E,5E,10E,15E]
        myYrange=[0,0]
     end
  endcase


  if keyword_set(FWHM) then begin
     GHW = sigma * 1.17741 ;; Gaussian Half Width at Half Max
     LHW = GHW * vdamp         ;; Lorentzian Half width at Half Max
     FWHMvoigt = (0.5346E * LHW + sqrt(0.2166E * LHW^2 + GHW^2)) * 2E
     sigmaString = 'FWHM = '+string(FWHMvoigt,format='(F6.1)')+' px '
  endif else begin
     sigmaString = cgGreek('sigma')+' = '+string(sigma,format='(F6.1)')+' px '
  endelse

  nlines = n_elements(sigma)
  farray = dblarr(nlines,npts)

  case 1 of 
     keyword_set(voigt): expr = 'voigt_slit(X,P[0],P[1],P[2])'
     keyword_set(avoigt): expr = 'vslit_approx(X,P[0],P[1],P[2])'
     else: expr = 'gauss_slit(X,P[0],P[1])'
  endcase

  for i=0l,nlines-1l do begin
     if keyword_set(differential) then begin
        farray[i,*] = expression_eval(expr,y,[H,sigma[i],vdamp]) - $
                      expression_eval(expr,y - starsep,[H,sigma[i],vdamp])
     endif else begin
        farray[i,*] = expression_eval(expr,y,[H,sigma[i],vdamp])
     endelse     
  endfor
  
  colorArray = [!p.color,mycol(['red','blue','dgreen'])]
  linestyleArray = [0,1,3,4]
  if keyword_set(differential) then begin
     myYtitle = 'Differential Transmission (%)'
  endif else myYtitle='Transmission (%)'

  myYextent = max(farray) - min(farray)
  if n_elements(myYRange) EQ 0 then begin
     myYrange = [min(farray),max(farray) + myYextent * 0.2] * 100E
  endif

  plot,y,farray[0,*] * 100E,ytitle=myYtitle,$
       xtitle='Position (pixels)',$
       ystyle=16+1,/nodata,yrange=myYrange
  for i=0l,nlines-1l do begin
     oplot,y,farray[i,*]*100E,color=colorArray[i],linestyle=lineStyleArray[i]
  endfor

  if keyword_set(psplot) then myCharsize=0.7 else myCharsize=1.0

  al_legend,sigmastring,color=colorArray,linestyle=linestyleArray,$
            charsize=myCharsize
  ;; second legend for explanatory text
  exptext = ['H = '+string(H,format='(F4.1)')+'px']
  
  if keyword_set(voigt) OR keyword_set(avoigt) then begin
     exptext = [exptext,'a = '+string(vdamp[0],format='(F6.2)')]
  endif
  if keyword_set(differential) then begin
     exptext = [exptext,cgGreek('Delta')+'y = '+string(starsep,format='(F6.2)')+'px']
  endif
  al_legend,exptext,/right,charsize=myCharsize

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

