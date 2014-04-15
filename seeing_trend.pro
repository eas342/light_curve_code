pro seeing_trend
;; This script plots seeing as a funcion of time, fits a polynomial
;; and then applies this to the time series

!p.multi=[0,1,2]

;; get the widths
restore,'data/prof_widths.sav'
y = widths[*,0]
nLength = n_elements(y)
x = findgen(nLength)
yerr = fltarr(nLength) + 0.1
expr = 'eval_poly(X,P)' ;; fitting function
startParams = float([1,fltarr(7)])
goodp = where(finite(y))
result = mpfitexpr(expr,x[goodp],y[goodp],yerr[goodp],startParams)
yfit = expression_eval(expr,x,result)

plot,x,y
oplot,x,yfit,color=mycol('red')
;; time serieso
restore,'data/specdata.sav'

normFl = transpose(binfl[0,*])
normFL = normFl / median(normFl)
yOrig = normFl
yfitNorm = yfit / median(yfit)
MixParam = -0.025D
yCorr = yOrig/((1D - MixParam) + MixParam * yfitNorm)
plot,x,yOrig,yrange=[0.97,1.02],psym=4
oplot,x,yCorr,color=mycol('red'),psym=5


!p.multi=0

end
