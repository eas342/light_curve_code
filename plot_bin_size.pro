pro plot_bin_size,psplot=psplot,scalephoton=scalephoton,$
                  wavlmode=wavlmode,photmode=photmode,$
                  custyrange=custyrange,custtitle=custtitle,$
                  nointerplots=nointerplots,tserRange=tserRange,$
                  overplotmode=overplotmode,photind=photind,$
                  showname=showname,photleg=photleg,smallleg=smallleg
;; plots the rms as a function of bin size
;; psplot -- makes eps, pdf and png plots
;; scalephoton -- scales the photon errors up to the measured RMS
;; wavlmode -- looks at binning as a function of wavelength, instead
;;             of time
;; photmode -- photometry mode only looks at one wavelength
;; custyrange -- allows you to set the Y range of the plot
;; custtitle -- a custom title for RMS vs wavelength plots
;; nointerplots - doesn't show the tim series plots. useful
;;                when trying to make multi-night plots of bin size vs time
;; overplotmode - overplot-only mode, does not re-create the plot
;; photind - if specified, the index for photometry, default is 0
;; showname - an optional name to show in overplot mode
;; photleg - if specified, show a simple legend for photometry
;; smallleg - shrink legend so it isn't so big in ps plot

if keyword_set(wavlmode) then begin
   bintarr = [5,7,9,12,15,25,30,60,100,0]
;   bintarr = [5,7,15]
endif else begin
   bintarr = [300,200,180,150,130,100,90,80,70,60,55,50,45,40,37,$
              35,33,30,27,25,23,20,18,15,14,13,12,10,9,8]
   ;; Skip time bins that are shorter than exposure times
   restore,'data/specdata.sav'
   maxpoints = n_elements(utgrid)/2
   allowedbins = where(bintarr LT maxpoints)
   bintarr = bintarr[allowedbins]
endelse

if n_elements(custtitle) EQ 0 then custtitle=''

ntbin = n_elements(bintarr)

if keyword_set(photmode) then begin
   selectwav = [0.5]
endif else selectwav = [1.26,1.43,2.14] ;; microns
nselectwav = n_elements(selectwav)
rmstbinfun = fltarr(ntbin,nselectwav)
bintsizeArr = fltarr(ntbin)
photonfun = fltarr(ntbin,nselectwav)
for i=0l,ntbin-1l do begin
   if keyword_set(wavlmode) then begin
      if bintarr[i] EQ 0 then begin
         compile_spec,/matchgrid,/readCurrent,/irafnoise
         restore,'data/specdata.sav'
         bintarr[i] = float(n_elements(lamgrid))
      endif else compile_spec,nwavbins=bintarr[i],/readCurrent,/irafnoise
      plot_tim_ser,noplots=nointerplots,custxrange=tserRange
   endif else plot_tim_ser,timebin=bintarr[i],noplots=nointerplots,custxrange=tserRange
   restore,'data/rmsdata.sav'

   ;; get the wavelength bin starts bingrid
   ;; get the wavelength bin middle bingridmiddle
   ;; get the wavelength bin sizes binsizes
   ;; the rms of off transit flux fracRMSarr
   ;; the time bin sizes tsizes
   ;; the photon noise fracPhotonarr
   if keyword_set(photmode) then begin
      if n_elements(photind) EQ 0 then photind=0
      wavInd = [photind]
   endif else begin
      tabinv,bingridmiddle,selectwav,wavInd
      wavInd = round(wavInd)
   endelse

   rmstbinfun[i,*] = fracRMSarr[wavInd]
   photonfun[i,*] = fracPhotonarr[wavInd]

   if not keyword_set(wavl) then begin
      bintsizeArr[i] = tsizes[0]
   endif
   if i mod 10 EQ 0 then print,i,' of ',ntbin,' done'
endfor

if keyword_set(wavlmode) then begin
   readcol,'data/wavelength_ranges.txt',StartWavArr,EndWavArr,$
           /silent,format='(F,F)',skipline=1
   myXtitle = 'Wavelength bin size (um)'
   bintdescrip = (EndWavArr[0] - StartWavArr[0]) / bintarr
endif else begin
   bintdescrip = (bintsizeArr * planetdat.period) * (24D * 60D)
   myXtitle = 'Time bin size (min)'
endelse 

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
   photonName = 'Minimum Noise x '+string(scaleFact,format='(F8.3)'+' ')
endif else begin
   photonName = 'Minimum Noise '
endelse
;myYrange=[min(photonfun) - max(rmstbinfun) * 0.4E,max(rmstbinfun) *
;                           1.2E]*100E
;myYrange=[2E-2,5E]
if n_elements(custyrange) EQ 0 then begin
   if keyword_set(wavlmode) then myYrange= [2E-2,5E] else myYrange=[1E-3,2E]
endif else myYrange=custyrange

if keyword_set(overplotmode) then begin
   showcol = mycol('red')
   if n_elements(showname) EQ 0 then showname=''
   oplot,bintdescrip,rmstbinfun[*,0]*100E,color=showcol
   xyouts,median(bintdescrip),median(rmstbinfun[*,0]*100E) * 2E,showname,color=showcol
endif else begin
   plot,bintdescrip,rmstbinfun[*,0]*100E,$
        xtitle=myXtitle,$
        ytitle='Out of Transit RMS (%)',$ ;/xlog,$
        title=custtitle,$
        yrange=myYrange,ystyle=1,/ylog
   showcol = !p.color
endelse

oplot,bintdescrip,photonfun[*,0]*100E,linestyle=2,color=showcol

if keyword_set(wavlmode) then begin
   wavbinnames = 'Near '+$
                 strtrim(string(bingrid[wavInd]+ binsizes[wavInd]/2E,format='(F8.2)'),1) +$
                 ' um'
endif else begin
   
   wavbinnames = string(bingrid[wavInd],format='(F7.2)') + 'um to'+$
                 string(bingrid[wavInd]+ binsizes[wavInd],format='(F7.2)')+$
                 'um'
endelse

if not keyword_set(photmode) then begin
   oplot,bintdescrip,rmstbinfun[*,1]*100E,color=mycol('red')
   oplot,bintdescrip,photonfun[*,1]*100E,color=mycol('red'),linestyle=2
   oplot,bintdescrip,rmstbinfun[*,2]*100E,color=mycol('blue')
   oplot,bintdescrip,photonfun[*,2]*100E,color=mycol('blue'),linestyle=2
   
   legend,'Std Dev '+wavbinnames,color=mycol(['black','red','blue']),linestyle=[0,0,0],/right
   legend,color=mycol(['black','red','blue']),replicate(photonName,3)+wavbinnames,$
          linestyle=replicate(2,3),/bottom,/left
endif else begin
   if keyword_set(photleg) then begin
      if keyword_set(smallleg) then begin
         legcharsize=0.5
      endif else legcharsize=1.0
      al_legend,['Measured Noise','Photon + read noise'],linestyle=[0,2],charsize=legcharsize
   endif
endelse

if keyword_set(psplot) then begin
   device, /close
   cgPS2PDF,plotprenm+'.eps'
   spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

end
