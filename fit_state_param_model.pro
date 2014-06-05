pro fit_state_param_model,psplot=psplot,fixspatial=fixspatial,$
                          psmooth=psmooth,widthfix=widthfix
;; Fits the state parameter model to the observed flux ratio time
;; series
;; psplot - saves a postscript plot
;; fixspatial - fixes the FWHM and voigt parameters in time (with the
;;              median value)
;; psmooth - smooth the parameters in time
;; widthfix - fixes the slit width at an estimated size

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
inputX = dblarr(7,npts)
;inputX[0,*] = statePstruct.position1 ;; Ack! Wrong position (spatial)
;inputX[1,*] = statePstruct.position2
inputX[0,*] = statePstruct.specshift1
inputX[1,*] = statePstruct.specshift2
inputX[2,*] = statePstruct.fwhm1
inputX[3,*] = statePstruct.fwhm2
inputX[4,*] = statePstruct.phase
inputX[5,*] = statePstruct.voigtdamp1
inputX[6,*] = statePstruct.voigtdamp2

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


;; Trim out the bad points, but make a copy of the original
originalX = inputX
originalY = y

AllGood = where(badArray EQ 0,ngood)
y = y[AllGood]
inputX = inputX[*,AllGood]
yerr = robust_sigma(y) + fltarr(ngood)

;; For now smooth the state parameters
if keyword_set(psmooth) then begin
   smoothP = [0,1,2,3,5,6]
   nsmoothP = n_elements(smoothP)
   for i=0l,nsmoothP-1l do begin
      inputX[smoothP[i],*] = smooth(inputX[smoothP[i],*],30)
   endfor
endif

if keyword_set(fixspatial) then begin
;inputX[1,*] = smooth(inputX[0,*],20)
;; For now set the FWHM & voigt Damp to be constant with time
   inputX[2,*] = dblarr(ngood) + median(inputX[2,*])
   inputX[3,*] = dblarr(ngood) + median(inputX[3,*])
   inputX[5,*] = dblarr(ngood) + median(inputX[5,*])
   inputX[6,*] = dblarr(ngood) + median(inputX[6,*])
end

Start = [0D,-1D,10.5D,1.0D,0D,0.3D,1.0D]
;Start = [1.5D,-1D,14D,0.985D,0D,0D,0D]
;fitexpr = 'gauss_slit(X[0,*] - P[0],P[2],X[2,*])/gauss_slit(X[1,*] - P[1],P[2],X[3,*]) *'+$
;          ' eval_legendre(X[4,*],P[3:4])'


nparams = n_elements(start)
pi = replicate({fixed:0, limited:[0,0], limits:[0.0E,0.0E]},nparams)
if keyword_set(widthfix) then pi[2].fixed = 1
pi[5].fixed = 1
;pi[5].limited = [1,1]
;pi[5].limits = [0E,1E]
fitexpr = 'vslit_approx(X[0,*] - P[0],P[2],X[2,*],X[5,*])/'+$
          'vslit_approx(X[1,*] - P[1],P[2],X[3,*]*P[6],X[6,*])'+$
          ' * eval_legendre(X[4,*],P[3:4])'
;fitexpr = 'voigt_slit(X[0,*] - P[0],P[2],X[2,*],0.4)/gauss_slit(X[1,*] - P[1],P[2],X[3,*]) *'+$
;          ' eval_legendre(X[4,*],P[3:4])'

result = mpfitexpr(fitexpr,inputX,y,yerr,start,parinfo=pi,perr=perr)

;; Show the model
;result = [1.5E,-1E,14E,0.985,-0.001D]
;result = [0E,-1E,10.5E,0.995D,-0.001D,0.3D]
ymodel = expression_eval(fitexpr,inputX,result)
xplot = statePstruct.phase

!p.multi = [0,1,2]
plot,xplot,y,yrange=[0.98,1.02],psym=4,$
     xtitle='Phase',ytitle='Normalized Flux'
oplot,xplot,ymodel,color=mycol('blue')

;; Print the parameters and errors
print,fitexpr
for i=0l,nparams-1l do begin
   print,'P['+strtrim(i,1)+'] = ',$
         strtrim(result[i],1),' +/- ',$
         strtrim(perr[i],1)
endfor

ycorrected = y / ymodel
print,'Robust sigma corrected = ',robust_sigma(y)
print,'Robust sigma corrected = ',robust_sigma(ycorrected)

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

yfullmodel = expression_eval(fitexpr,originalX,result)
yfullcorrected = originalY /yfullmodel

save,yfullcorrected,$
     filename='data/state_parameters/alt_tim_ser.sav'

end
