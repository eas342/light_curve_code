pro non_linearity_effect,psplot=psplot
;; Shows the corrections for Non-Linearity on the SpeX detector

  ;; set the plot
if keyword_set(psplot) then begin
   set_plot,'ps'
   !p.font=0
   plotprenm = 'plots/linearity/response_curve'
   device,encapsulated=1, /helvetica,$
          filename=plotprenm+'.eps'
   device,xsize=14, ysize=10,decomposed=1,/color
endif

restore,filepath('lc_coeff.sav', ROOT=file_dirname(file_which('SpeX.dat'),/MARK))

pedimage = lindgen(1,1024)
CorFac = lc_coeff[520,*,0]/mc_imgpoly(pedimage,lc_coeff[520,*,*])

plot,pedimage,CorFac,ystyle=16,$
     xtitle='Counts (DN)',$
     ytitle='Non-Linear Correction Factor'

x = transpose(pedimage)
y = transpose(CorFac)
;; Fit a line
;; dispose of Nans
goodp = where(finite(y))
x = x[goodp]
y = y[goodp]

rlinefit = robust_linefit(x,y,yfit)
oplot,x,yfit,color=mycol('red')

corHigh = rlinefit[0] + 600E * rlinefit[1]
corLow = rlinefit[0] + 588E * rlinefit[1]

print,'For a 2% transit depth at 600 DN'
print,' the correction factor goes from ',corHigh,' at DN =600'
print,' to ',CorLow,' at DN=588'


if keyword_set(psplot) then begin
   device, /close
   cgPS2PDF,plotprenm+'.eps'
   spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif


end
