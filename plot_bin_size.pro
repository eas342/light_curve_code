pro plot_bin_size,psplot=psplot,scalephoton=scalephoton
;; plots the rms as a function of bin size
;; psplot -- makes eps, pdf and png plots
;; scalephoton -- scales the photon errors up to the measured RMS


bintarr = [300,200,180,150,130,100,90,80,70,60,55,50,45,40,37,$
           35,33,30,27,25,23,20,18,15,14,13,12,10,9,8]
;bintarr = [200,6]
ntbin = n_elements(bintarr)

selectwav = [2.27,1.48,1.22] ;; microns
nselectwav = n_elements(selectwav)
rmstbinfun = fltarr(ntbin,nselectwav)
bintsizeArr = fltarr(ntbin)
photonfun = fltarr(ntbin,nselectwav)
for i=0l,ntbin-1l do begin
   plot_tim_ser,timebin=bintarr[i]
   restore,'data/rmsdata.sav'
   ;; get the wavelength bin starts bingrid
   ;; get the wavelength bin middle bingridmiddle
   ;; get the wavelength bin sizes binsizes
   ;; the rms of off transit flux fracRMSarr
   ;; the time bin sizes tsizes
   ;; the photon noise fracPhotonarr
   tabinv,bingridmiddle,selectwav,wavInd
   rmstbinfun[i,*] = fracRMSarr[round(wavInd)]
   photonfun[i,*] = fracPhotonarr[round(wavInd)]
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

if keyword_set(scalephoton) then begin
   scaleFact = 10E
   photonfun = photonfun * scaleFact
   photonName = 'Photon Noise x '+string(scaleFact,format='(F8.3)')
endif else begin
   photonName = 'Photon Noise'
endelse

plot,bintmin,rmstbinfun[*,0]*100E,$
     xtitle='Time bin size (min)',$
     ytitle='Out of Transit RMS (%)',$;/xlog,$
     yrang=[-0.2,0.4],ystyle=1
oplot,bintmin,photonfun[*,0]*100E,linestyle=2

wavbinnames = string(bingrid[wavInd],format='(F7.2)') + 'um to'+$
              string(bingrid[wavInd]+ binsizes[wavInd],format='(F7.2)')+$
              'um'

oplot,bintmin,rmstbinfun[*,1]*100E,color=mycol('red')
oplot,bintmin,photonfun[*,1]*100E,color=mycol('red'),linestyle=2
oplot,bintmin,rmstbinfun[*,2]*100E,color=mycol('blue')
oplot,bintmin,photonfun[*,2]*100E,color=mycol('blue'),linestyle=2

legend,wavbinnames+' Std Dev',color=mycol(['black','red','blue']),linestyle=[0,0,0],/right
legend,color=mycol(['black','red','blue']),replicate(photonName,3)+wavbinnames,$
       linestyle=replicate(2,3),/bottom,/left

if keyword_set(psplot) then begin
   device, /close
   cgPS2PDF,plotprenm+'.eps'
   spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

end
