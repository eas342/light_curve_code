pro show_training_results,pnum,psplot=psplot,custreject=custreject,$
                          domedian=domedian,custyrange=custyrange
;; Shows the fitted hyper-parameters from multi-train.pro
;; pnum -- the parameter index (0 for theta_0 1 for theta_1)
;; psplot -- save a postscript plot
;; custreject -- specifies the outlier rejection scheme to use when
;;               finding the mean hyperparameters
;; domedian -- find the median value of the best-fit parameters
;;             instead of the robust average
;; custyrange -- custom Y range

restore,'data/training_series.sav' ;; ParamArr and ParamTrue

;; set the plot
if keyword_set(psplot) then begin
   set_plot,'ps'
   !p.font=0
   plotprenm = 'plots/simulated_series/training_hist'
   device,encapsulated=1, /helvetica,$
          filename=plotprenm+'.eps'
   device,xsize=10, ysize=7,decomposed=1,/color
endif

yArray = ParamArr[pnum,*]
robustSig = robust_sigma(yArray)

yhist = histogram(yArray,binsize=RobustSig*0.3,locations=xhist)

paramsymbols = [cgGreek('theta')+'!D0!N',$
                cgGreek('theta')+'!D1!N',$
                cgGreek('sigma')]

plot,xhist,yhist,psym=10,$
     xtitle=paramsymbols[pnum],$
     ytitle='Counts',yrange=custyrange

;; Show the input value
oplot,ParamTrue[pnum]*[1,1],!y.crange,color=mycol('red'),linestyle=2

;; Show the average value & standard deviation
if n_elements(custreject) EQ 0 then custreject=3
meanVal = robust_mean(yArray,oreject=custreject,err=sigVal)
if keyword_set(domedian) then meanVal = median(yArray)
oplot,meanVal*[1,1],!y.crange,color=mycol('blue'),linestyle=0
oplot,(meanVal+sigVal)*[1,1],!y.crange,color=mycol('blue'),linestyle=1
oplot,(meanVal-sigVal)*[1,1],!y.crange,color=mycol('blue'),linestyle=1

al_legend,['True Input Value','Histogram','Robust Mean','Mean +/- Error'],$
       /right,linestyle=[2,0,0,1],color=[mycol('red'),!p.color,mycol(['blue','blue'])],$
       /clear


if keyword_set(psplot) then begin
   device, /close
   cgPS2PDF,plotprenm+'.eps'
   spawn,'convert -density 300% '+plotprenm+'.pdf '+plotprenm+'.png'
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

end
