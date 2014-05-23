pro fit_state_param_model,psplot=psplot
;; Fits the state parameter model to the observed flux ratio time
;; series
;; psplot - saves a postscript plot

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/spec_t_series/slit_model'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     device,xsize=15, ysize=10,decomposed=1,/color
     !p.thick=3.5
     !x.thick=2
     !y.thick=2
  endif
  
;; get the spectral info
restore,'data/specdata.sav'

;; get the spec list name for the observation
restore,'data/used_date.sav'

;; Read in the state parameters
restore,'data/state_parameters/full_parameters/'+$
   specfileListNamePrefix+'.sav'

npts = n_elements(statePstruct.phase)
inputX = dblarr(5,npts)
inputX[0,*] = statePstruct.position1
inputX[1,*] = statePstruct.position2
inputX[2,*] = statePstruct.fwhm1
inputX[3,*] = statePstruct.fwhm2
inputX[4,*] = statePstruct.phase

y = statePstruct.fluxratio
;; Normalize
y = y/median(y)

;; Reject all 10 sigma outliers
badArray = intarr(npts)
for i=0l,3+1l do begin
   if i EQ 4l then tempArray = y else tempArray = inputX[i,*]
   rsigma = robust_sigma(tempArray)
   goodp = where(abs(tempArray - median(tempArray)) LT 10E * rsigma,complement=badp)
   if badP NE [-1] then begin
      badArray[badp] = 1
   endif
endfor
;;; also try clipping ht e fist part of array
;earlyPt = where(statePstruct.phase LT 0.47)
;badArray[earlyPt] = 1

AllGood = where(badArray EQ 0)
y = y[AllGood]
inputX = inputX[*,AllGood]
Start = [0D,0D,10D,1.0D,0D]
;Start = [1.5D,-1D,14D,0.985D,0D,0D,0D]
fitexpr = 'gauss_slit(X[0,*] - P[0],P[2],X[2,*])/gauss_slit(X[1,*] - P[1],P[2],X[3,*]) *'+$
          ' eval_legendre(X[4,*],P[3:4])'
;fitexpr = 'voigt_slit(X[0,*] - P[0],P[2],X[2,*],0.4)/gauss_slit(X[1,*] - P[1],P[2],X[3,*]) *'+$
;          ' eval_legendre(X[4,*],P[3:4])'
result = mpfitexpr(fitexpr,inputX,y,fltarr(npts,4)+0.1,start)

;; Show the model
;result = [1.5E,-1E,14E,0.985,-0.001D]
;result = [1.5E,-1E,14E,0.985D,0D]
ymodel = expression_eval(fitexpr,inputX,result)
xplot = statePstruct.phase

!p.multi = [0,1,2]
plot,xplot,y,yrange=[0.98,1.02],psym=4,$
     xtitle='Phase',ytitle='Normalized Flux'

oplot,xplot,ymodel,color=mycol('blue')

ycorrected = y / ymodel

plot,xplot,ycorrected,yrange=[0.98,1.02],psym=4,$
     xtitle='Phase',ytitle='Corrected Flux'

  if keyword_set(psplot) then begin
     device, /close
;     cgPS2PDF,plotprenm+'.eps'
;     spawn,'convert -density 300% '+plotprenm+'.pdf'+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick=1
     !x.thick=1
     !y.thick=1
  endif
!p.multi=0

end
