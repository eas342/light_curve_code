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

x = binind[*,1,*]
y = binind[*,0,*]
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

linefit = robust_linefit(x,y,yfit)
oplot,x,yfit,color=mycol('lblue')
print,linefit[0],linefit[1]

resid = y - yfit
plot,x,resid,yrange=threshold(resid,mult=0),psym=4,$
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


!p.multi=0

end
