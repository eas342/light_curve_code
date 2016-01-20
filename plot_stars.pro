pro plot_stars,psplot=psplot,tryclean=tryclean,saveclean=saveclean,$
                  removelinear=removelinear,scalephoton=scalephoton,$
               flatten=flatten,smoothtemp=smoothtemp,choose1=choose1,$
               divide=divide,wavenum=wavenum,custXrange=custXrange,$
               showTelluric=showTelluric,custYrange=custYrange,normall=normall,$
               noNorm=noNorm,noLegend=noLegend,$
               showback=showback,directText=directText,custXmargin=custXmargin,$
               custYmargin=custYmargin,skipXtitle=skipXtitle,$
               choose2=choose2,digfilter=digfilter,biggerImage=biggerImage,$
               rmean=rmean,custtitle=custtitle,nobacknorm=nobacknorm
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
;; custX/Yrange -- allows you to input a custom X range instead of the
;;               full
;; normall -- normalize all spectra (both the host star and reference star)
;; noNorm - do not normalize spectra
;; showback -- show a background spectra
;; nobacknorm - don't normalize the background by itself so you
;;              can compare to the stars
;; directText -- show the text directly instead of with a legend
;; custX/Ymargin - used by double_spec to fine-tune margins
;; skipXtitle - skips X title and tick labels, for use by double specphot
;; digfilter -- apply a digital filter to the spectrum (convolve it)
;; noLegend - don't make a plot legend
;; biggerImage -- makes a bigger plot image (not so dense)
;; rmean - robust mean instead of median
;; custtitlte - custom plot title

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/individual_spectra/indspec'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     if keyword_set(biggerImage) then begin
        device,xsize=15, ysize=9,decomposed=1,/color
     endif else device,xsize=9, ysize=5,decomposed=1,/color
     !p.thick = 3.0
     !x.thick = 3.0
     !y.thick = 3.0
  endif


  ;; get the compiled spectroscopic data
  restore,'data/specdata.sav'

  restore,'data/used_date.sav'

  ;; get the time info
  restore,'data/timedata.sav'
  if tdataUseDate NE usedate then begin
     ;; Time data is not updated, so calculate it
     plot_tim_ser,/noplots
     restore,'data/timedata.sav'
  endif
  ;; divide all time series by the transit model

  tplot = find_phase(utgrid)
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
  backspec = fltarr(nwavs)
  ;; go through each wavelength and find the median
  cleanedcurve = fltarr(nwavs,ntime)
  cleanfactor = 4E
  oreject = 3E ;; outlier rejection for robust mean

  if n_elements(choose1) EQ 0 then begin
     for i=0l,nwavs-1l do begin
        if keyword_set(rmean) then begin
           hostspec[i] = robust_mean(flgrid[i,0,*],oreject=oreject)
           refspec[i] = robust_mean(flgrid[i,1,*],oreject=oreject)
           backspec[i] = robust_mean(backgrid[i,0,*],oreject=oreject)
        endif else begin
           hostspec[i] = median(flgrid[i,0,*])
           refspec[i] = median(flgrid[i,1,*])
           backspec[i] = median(backgrid[i,0,*])
        endelse
     endfor
  endif else begin
     hostspec = flgrid[*,0,choose1]
     refspec = flgrid[*,1,choose1]
     backspec = backgrid[*,0,choose1]
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

  multfac = 1E
;; Only show finite data, skip NANs
  goodp = where(finite(hostspec) EQ 1 and finite(refspec) EQ 1); and lamgrid LT 2.47)

  yhost = hostspec[goodp]
  yref = refspec[goodp] * multfac
  smoothsize = 20l
  if keyword_set(flatten) then begin
     yhost = yhost/smooth(yhost,smoothsize)
     yref = yref/smooth(yref,smoothsize)
  endif
  if keyword_set(digfilter) then begin
     yhost = convol(yhost,digital_filter(0.03,0.09,50,50))
;     yhost = convol(yhost,digital_filter(0.15,0.3,50,10))
     NormFac = max(yhost) - min(yhost)
     yhost = yhost/NormFac
  endif

  if keyword_set(flatten) then begin
     myYrange = [0.4,1.5]
     myYtitle='Normalized Flux + Offset'
  endif else begin
     myYrange=[0,1]
     myYtitle='Normalized Flux'
  endelse

  if n_elements(custXrange) EQ 0 then begin
     custXrange = [lamgrid[0],lamgrid[n_elements(lamgrid)-1l]]
     myXstyle=1
  endif else begin
     ;; Truncate to maximum possible
     custXrange = float(custXrange) ;; convert to float
     if custXrange[0] LT min(lamgrid) then custXrange[0] = min(lamgrid)
     if custxrange[1] GT max(lamgrid) then custXrange[1] = max(lamgrid)
     myXstyle=1
  endelse

  if n_elements(custYrange) NE 0 then myYrange = custYrange
  if n_elements(custXmargin) EQ 0 then custXmargin = [10,4]

;  middle80ref = threshold([yref,yhost],mult=0,high=0.95)
  wavthresh = threshold(lamgrid,mult=0,high=0.95)
  swavpt = where(lamgrid LT wavthresh[1])
  scaleData = [yref[swavpt],yhost[swavpt]]
  if keyword_set(nobacknorm) then begin
     scaleData = [scaleData,backspec[swavpt]]
  endif
     
  maxnorm = max(scaleData,/nan);middle80ref[1]

  case 1 of
     keyword_set(noNorm): ystarplot = yhost
     keyword_set(normall): begin
        middle80 = threshold(yhost,mult=0.05)
        ystarplot = yhost/middle80[1]
     end
     keyword_set(divide): begin
        ystarplot = yhost/yref
        yref = fltarr(n_elements(yref))
     end
     else: ystarplot = yhost/maxnorm
  endcase
  
  if keyword_set(skipXtitle) then begin
     myXtitle=''
     myXtickformat='(A1)'
  endif else begin
     myXtitle='Wavelength (um)'
     myXtickformat='(G8.2)'
  endelse

  plot,lamgrid[goodp],ystarplot,$
       xtitle=myXtitle,$
       ytitle=myYtitle,$
       yrange=myYrange,xrange=custXrange,xstyle=myXstyle,$
       xmargin=custXmargin,xtickformat=myXtickformat,$
       title=custtitle
;,xrange=[5600,6600],xstyle=1
;,ystyle=16;,xrange=[1.15,1.35]
;,xrange=[1.45,1.75]

  if keyword_set(normall) then begin
     middle60ref = threshold(yref,mult=0.05)
     yrefNorm = yref/middle60ref[1]

  endif else yrefNorm = yref/maxnorm
  oplot,lamgrid[goodp],yrefNorm,color=mycol('blue'),linestyle=3

  if n_elements(choose2) GT 0 then begin
     yhost2 = flgrid[*,0,choose2]/maxnorm
     yref2 = flgrid[*,1,choose2] /maxnorm
     oplot,lamgrid,yhost2,color=mycol('red')
     oplot,lamgrid,yref2,color=mycol('dgreen'),linestyle=3
  endif

;  Show the smoothed spec
;  oplot,lamgrid[goodp],smooth(yhost,smoothsize),color=mycol('red')

  ;; Show error
  if n_elements(choose1) GT 0 then begin
     oplot,lamgrid[goodp],yerr[goodp]
  endif

  reftext = 'Ref Star'
  if not keyword_set(flatten) then reftext = reftext + ' * '+string(multfac,format='(F8.4)')

  case 1 of 
     keyword_set(showtelluric): name3 = 'Telluric Model'
     keyword_set(showback): name3 = 'Background'
     else: name3 = 'G0 V Template'
  endcase

;  legend,['Host Star',reftext,name3],$
;         color=mycol(['black','blue','red']),/right,linestyle=[0,3,4]
  if keyword_set(directText) then begin
     xsize = !x.crange[1] - !x.crange[0]
     ysize = !y.crange[1] - !y.crange[0]
     xrefText = 0.65E * xsize+!x.crange[0]
     yrefText = 0.7E * ysize + !y.crange[0]
     xyouts,0.4E * xsize+!x.crange[0],0.85E * ysize + !y.crange[0],$
            'Planet Host',color=!p.color,charsize=0.7
     xyouts,xrefText,yrefText,$
            'Reference Star',color=mycol('blue'),charsize=0.7
     xyouts,0.07E * xsize+!x.crange[0],0.2E * ysize + !y.crange[0],$
            name3,color=mycol('red'),charsize=0.7
     midpt = n_elements(lamgrid)/2
     dataperYpix = (!y.crange[1] - !y.crange[0]) /((!y.window[1] - !y.window[0]) * !d.y_vsize)
     ybump = !D.Y_CH_SIZE * 0.5E * dataperYpix * 0.7E
     oplot,[lamgrid[goodp[midpt]],xrefText],[yrefNorm[midpt],yrefText + ybump],color=mycol('blue')
  endif else begin
     case 1 of 
        keyword_set(noLegend): print,''
        keyword_set(choose2): begin
           legNames = ['Planet Host (Img '+strtrim(choose1,1)+')',$
                       'Planet Host (Img '+strtrim(choose2,1)+')',$
                       'Ref Star (Img '+strtrim(choose1,1)+')',$
                       'Ref Star (Img '+strtrim(choose2,1)+')']
           al_legend,legNames,$
                  color=[!p.color,mycol(['red','blue','dgreen'])],/right,$
                  linestyle=[0,0,3,3],charsize=0.45
        end
        else: begin
           al_legend,['Planet Host','Reference Star',name3],$
                  color=[!p.color,mycol(['blue','red'])],/right,linestyle=[0,3,4],$
                  charsize=0.7
        end
     endcase
  endelse

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

  if not keyword_set(showtelluric) then begin
     oplot,tempwavl[goodp],ytemp,color=mycol('red'),linestyle=4 ;* tempwavl^(1/2)
  endif
  


;  stop

;  Show the smoothed spec
;  polyfit = poly_fit(tempwavl[goodp],ytemp,8,yfit=ytemppoly)
;  oplot,tempwavl[goodp],smooth(ytemp,tempsmooth),color=mycol('green')
;  oplot,tempwavl[goodp],ytemppoly,color=mycol('green')

  if keyword_set(showTelluric) then begin
     restore,'data/telluric/mauna_kea_trans_gemini.sav'
     ;; Interpolate the 
     oplot,wavlt,smooth(trans,2000)/max(smooth(trans,2000)),color=mycol('red')
  endif

  if keyword_set(showback) then begin
     if keyword_set(nobacknorm) then bnormFac = maxnorm else bnormFac=(10E * median(backspec))
     oplot,lamgrid,backspec/bnormFac,color=mycol('red')

  endif

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick = 1.0
     !x.thick = 1.0
     !y.thick = 1.0
 endif

end
