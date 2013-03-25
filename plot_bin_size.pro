pro plot_bin_size,psplot=psplot
;; plots the rms as a function of bin size
;; psplot -- makes eps, pdf and png plots


bintarr = [200,150,100,80,60,50,40,35,30,25,20,15,13,10,8]
;bintarr = [200,6]
ntbin = n_elements(bintarr)

selectwav = [2.27,1.22] ;; microns
nselectwav = n_elements(selectwav)
rmstbinfun = fltarr(ntbin,nselectwav)
bintsizeArr = fltarr(ntbin)
for i=0l,ntbin-1l do begin
   plot_tim_ser,timebin=bintarr[i]
   restore,'data/rmsdata.sav'
   ;; get the wavelength bin starts bingrid
   ;; get the wavelength bin middle bingridmiddle
   ;; get the wavelength bin sizes binsizes
   ;; the rms of off transit flux fracRMSarr
   ;; the time bin sizes tsizes
   tabinv,bingridmiddle,selectwav,wavInd
   rmstbinfun[i,*] = fracRMSarr[round(wavInd)]
   bintsizeArr[i] = tsizes[0]
endfor

bintmin = (bintsizeArr * planetdat.period) * (24D * 60D)

;; set the plot
if keyword_set(psplot) then begin
   set_plot,'ps'
   !p.font=0
   plotprenm = 'plots/binsize/tbin_size'
   device,encapsulated=1, /helvetica,$
          filename=plotprenm+'.eps'
   device,xsize=14, ysize=10,decomposed=1,/color
endif



plot,bintmin,rmstbinfun[*,0]*100E,$
     xtitle='Time bin size (min)',$
     ytitle='Off Transit, Linearly-Detrended RMS (%)'

wavbinnames = string(bingrid[wavInd],format='(F10.3)') + 'um to'+$
              string(bingrid[wavInd]+ binsizes[wavInd],format='(F10.3)')+'um'

oplot,bintmin,rmstbinfun[*,1]*100E,color=mycol('red')


legend,wavbinnames,color=mycol(['black','red']),linestyle=[0,0]

if keyword_set(psplot) then begin
   device, /close
   cgPS2PDF,plotprenm+'.eps'
   spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

end
