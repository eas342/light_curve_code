pro fit_state_param_model,psplot=psplot,fixspatial=fixspatial,$
                          psmooth=psmooth,widthfix=widthfix,$
                          avgtwo=avgtwo,custyrange=custyrange,$
                          transit=transit,fixpos=fixpos,$
                          measuredDiff=measuredDiff,slitPos=slitPos,$
                          secondary=secondary,determslit=determslit
;; Fits the state parameter model to the observed flux ratio time
;; series
;; psplot - saves a postscript plot
;; fixspatial - fixes the FWHM and voigt parameters in time (with the
;;              median value)
;; psmooth - smooth the parameters in time
;; widthfix - fixes the slit width at an estimated size
;; avgtwo - Average the two stars' parameters
;; custyrange - set the y range of the plots
;; transit - fit light curve with transit
;; secondary - fit light curve w/ secondary eclipse
;; measuredDiff - use this to set the stars' separation at the
;;                measured value from the cross-correlation
;; slitPos - use a given fixed slit position instead of fitting for it
;; determslit -- use a determinate slit model (no fitted parameters)

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

SlitH = 10.5E
if n_elements(custyrange) EQ 0 then custyrange=[0.98,1.02]
  
;; get the spectral info
restore,'data/specdata.sav'

;; get the spec list name for the observation
restore,'data/used_date.sav'

;; get the planet info
readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
        skipline=1
planetdat = create_struct('null','')
for l=0l,n_elements(info)-1l do begin
   planetdat = create_struct(planetdat,info[l],data[l])
endfor


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

if keyword_set(avgtwo) then begin
   star1ind = [2,5]
   star2ind = [3,6]
   nstarts = n_elements(star1ind)
   for i=0l,nstarts-1l do begin
      avgX = (inputX[star1ind[i],*] + inputX[star2ind[i],*])/2E
      inputX[star1ind[i],*] = avgX
      inputX[star2ind[i],*] = avgX
   endfor
endif

y = statePstruct.fluxratio
;; Normalize
y = y/median(y)

;; Reject all n sigma outliers in state parameters
sigthreshS=5E
badArray = intarr(npts)
for i=0l,6l do begin
   if i EQ 4l then tempArray = y else tempArray = inputX[i,*]
   rsigma = robust_sigma(tempArray)
   goodp = where(abs(tempArray - median(tempArray)) LT sigthreshS * rsigma,complement=badp)
   if badP NE [-1] then begin
      badArray[badp] = 1
   endif
endfor

;;clip all n-sigma outliers in flux
sigthreshF=5E
rsigmaY = robust_sigma(y)
goodp = where(abs(y - median(y)) LT sigthreshF * rsigmaY,complement=badp)
badArray[badp] = 1

;;; also try clipping ht e fist part of array
;earlyPt = where(statePstruct.phase LT 0.47)
;badArray[earlyPt] = 1


;; Trim out the bad points, but make a copy of the original
originalX = inputX
originalY = y
noriginal = n_elements(originalY)

AllGood = where(badArray EQ 0,ngood,complement=allbad)
y = y[AllGood]
inputX = inputX[*,AllGood]
yerr = robust_sigma(y) + fltarr(ngood)

;; For now smooth the state parameters
if keyword_set(psmooth) then begin
   smoothP = [0,1,2,3,5,6]
   nsmoothP = n_elements(smoothP)
   if psmooth EQ 1 then smoothsize = 30 else smoothsize=psmooth
   for i=0l,nsmoothP-1l do begin
      inputX[smoothP[i],*] = smooth(inputX[smoothP[i],*],smoothsize)
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

if keyword_set(fixpos) then begin
   inputX[0,*] = 0E;dblarr(ngood) + median(inputX[0,*])
   inputX[1,*] = 0E;dblarr(
endif

Start = [0D,0D,slitH,1.0D,0D,0.0D,1.0D]
;Start = [1.5D,-1D,14D,0.985D,0D,0D,0D]
;fitexpr = 'gauss_slit(X[0,*] - P[0],P[2],X[2,*])/gauss_slit(X[1,*] - P[1],P[2],X[3,*]) *'+$
;          ' eval_legendre(X[4,*],P[3:4])'

if keyword_set(transit) then begin
   u1parm = 0.1E
   u2parm = 0.0E
   ;; Include the primary transit in model
   start = [start,0E,planetdat.p,planetdat.b_impact,u1parm,u2parm,$
            planetdat.a_o_rstar]
endif
if keyword_set(secondary) then begin
   u1parm = 0E
   u2parm = 0E
   start = [start,0E,0.001E,planetdat.b_impact,u1parm,u2parm,$
           planetdat.a_o_rstar]
endif

nparams = n_elements(start)
pi = replicate({fixed:0, limited:[0,0], limits:[0.0E,0.0E],tied:''},nparams)
if keyword_set(widthfix) then pi[2].fixed = 1
pi[5].fixed = 1

if keyword_set(transit) or keyword_set(secondary) then begin
   pi[9].fixed=1 ;; b_impact
   pi[11].fixed=1 ;; quadratic limb darkening
   pi[12].fixed=1 ;; a/R*
   if keyword_set(secondary) then pi[10].fixed=1 ;; linear limb darkening
endif

if keyword_set(measuredDiff) then begin
   pi[1].tied = 'P[0]'
endif

if keyword_set(determslit) then begin
   pi[0].fixed=1 ;; slit pos
   pi[1].fixed=1 ;; slit pos 2
   pi[2].fixed=1 ;; slit width
;   pi[4].fixed=1 ;; let's start with no linear
;   pi[6].fixed=1 ;; no offset in stellar profile
endif

if n_elements(slitPos) NE 0 then begin
   start[0] = slitPos
   pi[0].fixed = 1
endif

;; Make ethe stars be within the slit
positionInd = [0,1]
pi[positionInd].limited=[1,1]
pi[positionInd].limits = [-slitH,slitH]

;pi[5].limited = [1,1]
;pi[5].limits = [0E,1E]
fitexpr = 'vslit_approx(X[0,*] - P[0],P[2],X[2,*],X[5,*])/'+$
          'vslit_approx(X[1,*] - P[0] - P[1],P[2],X[3,*]*P[6],X[6,*])'+$
          ' * eval_legendre(X[4,*],P[3:4])'
;fitexpr = 'voigt_slit(X[0,*] - P[0],P[2],X[2,*],0.4)/gauss_slit(X[1,*] - P[1],P[2],X[3,*]) *'+$
;          ' eval_legendre(X[4,*],P[3:4])'

if keyword_set(transit) then begin
   fitexpr = fitexpr + ' * quadlc(X[4,*] - P[7],P[8],P[9],P[10],P[11],P[12])'
endif

if keyword_set(secondary) then begin
   fitexpr = fitexpr + ' * sec_eclipse(X[4,*] - P[7],P[8],P[9],P[12])'
endif

result = mpfitexpr(fitexpr,inputX,y,yerr,start,parinfo=pi,perr=perr)

;; Show the model
;result = [1.5E,-1E,14E,0.985,-0.001D]
;result = [0E,-1E,10.5E,0.995D,-0.001D,0.3D]
if keyword_set(secondary) then begin
   tempresult = result
;   tempresult[8] = 0E ;; sometimes I like to see the effect of the transit
   ymodel = expression_eval(fitexpr,inputX,tempresult)
endif else begin
   ymodel = expression_eval(fitexpr,inputX,result)
endelse
xplot = statePstruct.phase

!p.multi = [0,1,2]
plot,xplot,y,yrange=custyrange,psym=4,$
     xtitle='Phase',ytitle='Normalized Flux'
if n_elements(allbad) NE [-1] then begin
   oplot,originalX[allbad],originalY[allbad],color=mycol('red'),psym=2
endif
oplot,xplot,ymodel,color=mycol('blue')

;; Print the parameters and errors
print,fitexpr
for i=0l,nparams-1l do begin
   print,'P['+strtrim(i,1)+'] = ',$
         strtrim(result[i],1),' +/- ',$
         strtrim(perr[i],1)
endfor

ycorrected = y / ymodel
yresid = y - ymodel
print,'Robust sigma corrected = ',robust_sigma(y)
print,'Robust sigma corrected = ',robust_sigma(ycorrected)

plot,xplot,yresid,yrange=custyrange - 1.0E,psym=4,$
     xtitle='Phase',ytitle='Residual'

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

;yfullmodel = expression_eval(fitexpr,originalX,result)
;; Since we smoothed, averaged and modified inputX, I think it makes
;; sense to use the corrected time series and then !values.f_nan for
;; the points not included
yfullcorrected = !values.f_nan + fltarr(noriginal)
yfullcorrected[allgood] = ycorrected

save,yfullcorrected,$
     filename='data/state_parameters/alt_tim_ser.sav'
;stop
end
