pro two_star_corr,renorm=renorm,psplot=psplot,labelwav=labelwav
;; Compares the flux ratio of the two stars (assuming no transit)
;; Maybe it's not a 1:1 ratio. Why didn't I start with this!!??
;; renorm  - re-normalize the stars' fluxes by the
;;           target/reference ratio spectrum so they should match up
;; psplot - postscript plot
;; labelwav - label the wavelengths

if keyword_set(psplot) then begin
   set_plot,'ps'
   !p.font=0
   plotprenm = 'plots/linearity/two_star_corr'
   device,encapsulated=1, /helvetica,$
          filename=plotprenm+'.eps'
   device,xsize=20, ysize=15,decomposed=1,/color
   !p.thick=3.5
   !x.thick=2
   !y.thick=2
endif

restore,'data/specdata.sav'

medianSpec = median(binfl,dimension=2)
;; Re-normalize by the median spectrum to take out difference in color
;; between the two stars
if keyword_set(renorm) then begin
   for i=0l,nwavbins-1l do begin
;      binind[i,0,*] = binind[i,0,*] / median(binind[i,0,*]) * median(binind[i,1,*])
      binind[i,1,*] = binind[i,1,*] * medianSpec[i]
   endfor
endif

x2d = binind[*,1,*] ;; reference star individual flux
y2d = binind[*,0,*] ;; target star individual flux
;; make 1 d arrays
x = reform(x2d,n_elements(x2d))
y = reform(y2d,n_elements(y2d))
ntime = n_elements(utgrid)

!p.multi=[0,1,3]

xGlobalRange=threshold(x,mult=0)

plot,x,y,ystyle=16,psym=4,$
     xtitle='Ref Flux',$
     ytitle='Host Flux',$
     xrange=xglobalRange,yrange=threshold(y,mult=0.1),$
     charsize=2,/nodata

colorChoices = myarraycol(nwavbins)

for i=0l,nwavbins-1l do begin
   oplot,x2d[i,0,*],y2d[i,0,*],color=colorChoices[i],psym=4
endfor

if keyword_set(labelwav) then begin
   al_legend,[string(bingrid[0],format='(F5.2)')+' to '+$
             string(bingrid[0] + binsizes[0],format='(F5.2)')+' '+cgGreek('mu')+'m']
endif
;linefit = robust_linefit(x,y,yfit)
;oplot,x,yfit,color=mycol('lblue')
;print,linefit[0],linefit[1]

nIter = 3
expr = 'eval_poly(X,P)' ;; fitting function
yerr = fltarr(n_elements(y)) + robust_sigma(y)
sigmaThresh = 5E
for i=0l,nIter-1l do begin
   startParams = [0E,1E]
   nparams = n_elements(startParams)
   pi = replicate({fixed:0, limited:[0,0], limits:[0.0E,0.0E],tied:''},nparams)
   pi[0].fixed = 1
   if i EQ 0 then begin
      goodp = where(abs(y - median(y)) LE robust_sigma(y) * sigmaThresh,complement=badP)
   endif else begin
      goodp = where(abs(y - yfit) LE robust_sigma(y - yfit) * sigmaThresh,complement=badP)
   endelse
   result = mpfitexpr(expr,x[goodp],y[goodp],yerr[goodp],startParams,parinfo=pi)
;   yfit = expression_eval(expr,sortx,result)
;   oplot,sortx,yfit,color=mycol('lblue')
   yfit = expression_eval(expr,x,result)

endfor
sortx = x[sort(x)]
oplot,sortx,expression_eval(expr,sortx,result),color=mycol('lblue')
if badP NE [-1] then oplot,x[badP],y[badp],psym=4,color=mycol('red')

resid = y - yfit
plot,x,resid,yrange=threshold(resid,mult=0.3,low=0.1,high=0.9),psym=4,$
     xtitle='Ref Flux',$
     ytitle='Residual',$
     charsize=2,xrange=xGlobalRange

plot,x,y/x,psym=4,$
     ystyle=16,$
     xtitle='Ref Flux',$
     ytitle='Host/Ref Ratio',$
     yrange=threshold(y/x,mult=0),$
     xrange=xGlobalRange,$
     charsize=2

altY = y/yfit ;; polynomial correction
PosPt = where(x GT 0);; positive points
regularY = fltarr(n_elements(y)) + !values.f_nan
regularY[PosPt] = y[PosPt] / x[PosPt]

save,altY,regularY,$
     filename='data/alt_tim_ser.sav'

!p.multi=0

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
