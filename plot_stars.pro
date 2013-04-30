pro plot_stars,psplot=psplot,tryclean=tryclean,saveclean=saveclean,$
                  removelinear=removelinear,scalephoton=scalephoton,$
               flatten=flatten,smoothtemp=smoothtemp,choose1=choose1,$
               divide=divide,wavenum=wavenum,custXrange=custXrange
;; Plots the reference star and planet host
;; spectrum
;; psplot -- makes a postscript plot of the RMS spectrum
;; tryclean - tries to clean up the spectrum on a wavelength by
;;            wavelength basis
;; saveclean - saves the white light curve for analysis by plot_tim_ser
;; removelinear -- fit each curve to a line first before finding the RMS
;; scalephoton=scalephoton
;; flatten -- flattens the spectrum by dividing by a smoothed version
;; smoothtemp -- smooths the template by smoothtemp points
;; choose1 -- plots a single spectrum (instead of the median)
;; divide -- divides the stellar host by reference star
;; wavnum -- converts wavelength to wave number (cm^-1)
;; custXrange -- allows you to input a custom X range instead of the
;;               full

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/individual_spectra/indspec'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif


  ;; get the compiled spectroscopic data
  restore,'data/specdata.sav'

  ;; get the time info
  restore,'data/timedata.sav'

  ;; divide all time series by the transit model
  ymodel = quadlc(tplot,planetdat.p,planetdat.b_impact,$
                  u1parm,u2parm,planetdat.a_o_rstar)
  divsize = size(DivSpec)
  replicatedmodel = rebin(ymodel,divsize[3],divsize[2],divsize[1])
  rebinmodel = transpose(replicatedmodel,[2,1,0])
  divbycurve = flgrid[*,0,*] / rebinmodel

  if keyword_set(wavenum) then lamgrid = 1E4/lamgrid

  nwavs = n_elements(lamgrid)
  ntime = n_elements(tplot)
  hostspec = fltarr(nwavs)
  refspec = fltarr(nwavs)
  ;; go through each wavelength and find the median
  cleanedcurve = fltarr(nwavs,ntime)
  cleanfactor = 4E
  if n_elements(choose1) EQ 0 then begin
     for i=0l,nwavs-1l do begin
;        hostspec[i] = median(flgrid[i,0,*])
;        refspec[i] = median(flgrid[i,1,*])
        hostspec[i] = mean(flgrid[i,0,*],/nan)
        refspec[i] = mean(flgrid[i,1,*],/nan)
;     if i EQ 50l then stop
     endfor
  endif else begin
     hostspec = flgrid[*,0,choose1]
     refspec = flgrid[*,1,choose1]
     yerr = ErrGrid[*,0,choose1]
  endelse
  

  if keyword_set(tryclean) or keyword_set(removelinear) then begin
     for i=0l,nwavs-1l do begin
        if keyword_set(removelinear) then cleanedcurve[i,*] = divbycurve[i,0,*] else begin
           cleanedcurve[i,*] = DivSpec[i,0,*]
        endelse
        badp = where(abs(divbycurve[i,0,*] - median(divbycurve[i,0,*])) GT sigarray[i] * divbycurve[i,0,*] * 6E)
        if badp NE [-1] then cleanedcurve[i,badp] = !values.f_nan
     endfor
  endif

  if keyword_set(removelinear) then begin
     for i=0l,nwavs-1l do begin
        ;; divide by the line to flatten out
        goodp = where(finite(cleanedcurve[i,*]) EQ 1,ngoodp)
        if ngoodp GT 10 then begin
           rlinefit = robust_linefit(tplot[goodp],cleanedcurve[i,goodp],yfit)
           yflat = cleanedcurve[i,goodp] / yfit
           ;; You don't have to divide by the median because
           ;; the division by the line fit already does that
           sigarray[i] = robust_sigma(yflat)
        endif else sigarray[i] = !values.f_nan
        
     endfor
  endif

  if keyword_set(scalephoton) then begin
     scalefact = 4.2E
     divspecE[*,0,0] = divspecE[*,0,0] * scalefact
     photname = 'Photon Error x '+string(scalefact,format='(F8.2)')
  endif else photname = 'Photon Error'

  multfac = 0.42E
;; Only show finite data, skip NANs
  goodp = where(finite(hostspec) EQ 1 and finite(refspec) EQ 1)

  yhost = hostspec[goodp]
  yref = refspec[goodp] * multfac
  smoothsize = 20l
  if keyword_set(flatten) then begin
     yhost = yhost/smooth(yhost,smoothsize)
     yref = yref/smooth(yref,smoothsize)
  endif

  if keyword_set(flatten) then begin
     myYrange = [0.4,1.5]
  endif else myYrange=[0,max(yhost)*1.3]

  if n_elements(custXrange) EQ 0 then custXrange = [0,0]

  plot,lamgrid[goodp],yhost,$
       xtitle='Wavelength (um)',$
       ytitle='Flux (e!E-!N)',$
       yrange=myYrange,xrange=custXrange
;,xrange=[5600,6600],xstyle=1
;,ystyle=16;,xrange=[1.15,1.35]
;,xrange=[1.45,1.75]
  oplot,lamgrid[goodp],yref - 0.2,color=mycol('blue'),linestyle=3
  
;  Show the smoothed spec
;  oplot,lamgrid[goodp],smooth(yhost,smoothsize),color=mycol('red')

  ;; Show error
  if n_elements(choose1) GT 0 then begin
     oplot,lamgrid[goodp],yerr[goodp]
  endif

  reftext = 'Ref Star'
  if not keyword_set(flatten) then reftext = reftext + ' * '+string(multfac,format='(F8.4)')
  legend,['Host Star',reftext,'G0 V Template'],$
         color=mycol(['black','blue','red']),/right,linestyle=[0,3,4]

  ;; Read in a stellar template and over-plot
  readcol,'star_templates/g0v_template_irtf.txt',tempwavl,tempfl,tempflerr,$
          skipline=47,format='(F,F,F)'

  if keyword_set(wavenum) then tempwavl = 1E4/tempwavl
;  multfac2 = 3.5E15
  multfac2 = 2.5E15

  goodp = where(tempfl NE -999) ;; the NaNs are treated as -999

  ytemp = tempfl[goodp] * multfac2
  tempsmooth = 200l
  if keyword_set(flatten) then begin
     polyfit = poly_fit(tempwavl[goodp],ytemp,8,yfit=ytemppoly)
     ytemp = ytemp/ytemppoly - 0.4
  endif
  
  if n_elements(smoothtemp) NE 0 then ytemp = smooth(ytemp,smoothtemp)

  oplot,tempwavl[goodp],ytemp,color=mycol('red'),linestyle=4 ;* tempwavl^(1/2)

  


;  stop

;  Show the smoothed spec
;  polyfit = poly_fit(tempwavl[goodp],ytemp,8,yfit=ytemppoly)
;  oplot,tempwavl[goodp],smooth(ytemp,tempsmooth),color=mycol('green')
;  oplot,tempwavl[goodp],ytemppoly,color=mycol('green')


  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

  if keyword_set(tryclean) then begin
;     wait,1
     
     ;; Make a special broadband light curve
     specialpt = where(sigarray LT 0.01,nspecial)

     ;; Weight the points by the rms
;     weights=1/sigarray^2
     weights=1/sigarray^2
     ;; discount wavelengths and/or pixels that vary wildly
;     badp = where(sigarray GT 0.01 OR lamgrid LT 0.87 OR lamgrid GT 2.4)
;     badp = where(sigarray GT 0.20)
     if badp NE [-1] then weights[badp] = 0.0E
     weightsCopy = rebin(weights,nwavs,ntime)
     ;; Set all non-finite values of the cleaned curve to have NANS in the weights
     nonfinite = where(finite(cleanedcurve) EQ 0)
     weightsCopy[nonfinite] = !values.f_nan
     ;; weighted sum normalized by the weights for that spectrum
     ;; this accounts for missing data (that has NANs)
     combinedpt = total(cleanedcurve * weightsCopy,1,/nan)/total(weightsCopy,1,/nan)
     
     

     combinedpt2 = total(cleanedcurve[specialpt,*],1)/float(nspecial)
     plot,tplot,combinedpt2,ystyle=16,psym=4
;     plot,tplot,combinedpt2,ystyle=16,psym=4,yrange=[0.40,0.45]
;     oplot,tplot,combinedpt+0.01,psym=5,color=mycol('green')
;     plot,tplot,combinedpt,psym=5,color=mycol('green'),ystyle=16


     ;try out straight avg over wavelength range
     goodrange = where(lamgrid GE 0.9 and lamgrid LT 2.4,ngoodrange)
     combinedpt2 = total(cleanedcurve[goodrange,*],1)/float(ngoodrange)

;     plot,tplot,combinedpt2,ystyle=16,psym=4
;
     if keyword_set(saveclean) then begin
        
     endif
  endif



end
