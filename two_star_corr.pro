pro two_star_corr,renorm=renorm
;; Compares the flux ratio of the two stars (assuming no transit)
;; Maybe it's not a 1:1 ratio. Why didn't I start with this!!??

restore,'data/specdata.sav'

;; Re-normalize so that the median host flux for the wavelength bin
;; equals the median reference star flux
;; this is kind of pointless so I should take it out
if keyword_set(renorm) then begin
   for i=0l,nwavbins-1l do begin
      binind[i,0,*] = binind[i,0,*] / median(binind[i,0,*]) * median(binind[i,1,*])
   endfor
endif

x = binind[*,1,*] ;; reference star individual flux
y = binind[*,0,*] ;; target star individual flux
;; make 1 d arrays
x = reform(x,n_elements(x))
y = reform(y,n_elements(y))


!p.multi=[0,1,3]

xGlobalRange=threshold(x,mult=0)

plot,x,y,ystyle=16,psym=4,$
     xtitle='Ref Flux',$
     ytitle='Host Flux',$
     xrange=xglobalRange,yrange=threshold(y,mult=0.1),$
     charsize=2

;linefit = robust_linefit(x,y,yfit)
;oplot,x,yfit,color=mycol('lblue')
;print,linefit[0],linefit[1]

nIter = 3
expr = 'eval_poly(X,P)' ;; fitting function
yerr = fltarr(n_elements(y)) + robust_sigma(y)
sigmaThresh = 5E
for i=0l,nIter-1l do begin
   startParams = [1,0,0]
   if i EQ 0 then begin
      goodp = where(abs(y - median(y)) LE robust_sigma(y) * sigmaThresh,complement=badP)
   endif else begin
      goodp = where(abs(y - yfit) LE robust_sigma(y - yfit) * sigmaThresh,complement=badP)
   endelse
   result = mpfitexpr(expr,x[goodp],y[goodp],yerr[goodp],startParams)
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

end
