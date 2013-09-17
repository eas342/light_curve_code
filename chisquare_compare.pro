pro chisquare_compare
;; Calculates the chi-squared of spectral models with the data stored
;; in avg rp_rs

readcol,'radius_vs_wavelength/avg_rp_rs.txt',wavl,wavlsize,rad,rade,$
        skipline=1,format='(F,F,F,F)'

;; Add the Bean 2009 result
CorWav = 0.65         ;; microns, approximately
CorWid = 0.20         ;; microns, approx
CorRad = 0.1433 ;;Rp/R*
CorErr = 0.0010
wavl = [CorWav,wavl]
wavlsize = [CorWid,wavlsize]
rad = [CorRad,rad]
rade = [CorErr,rade]

wavlwidth = wavlsize

readcol,'../models/fortney_g10mps_2500K_isothermal.csv',theowav,theorad,$
        skipline=6,format='(F,F)'


ntheo=n_elements(theorad)
binModel = avg_series(theowav,theorad,fltarr(ntheo)+0.2E,wavl-wavlwidth/2E,wavlwidth,weighted=0)
save,wavl,binModel,filename='data/binned_model.sav'

nwavs = n_elements(wavl)

start1 = [0.1D]
expr1 = 'model_evaluate(X,P[0])'
result = mpfitexpr(expr1,wavl,rad,rade,start1)
ymod1 = expression_eval(expr1,wavl,result)
dof1 = nwavs - n_elements(start1)
chisq1 = total(((rad - ymod1)/rade)^2)/dof1
print,'Literature Model Chisq/DOF= ',chisq1

start2 = [0.1D]
expr2 = 'P[0]'
resultFlat = mpfitexpr(expr2,wavl,rad,rade,start2)
ymod2 = replicate(resultFlat,nwavs)
dof2 = nwavs - n_elements(start2)
chisq2 = total(((rad - ymod2)/rade)^2)/dof2
print,'Flat Model ChisQ/DOF = ',chisq2

plot,wavl,rad,/nodata,ystyle=16
oploterror,wavl,rad,wavlwidth/2E,rade,psym=3
oplot,wavl,ymod1,color=mycol('yellow')
oplot,wavl,ymod2,color=mycol('lblue')



end
