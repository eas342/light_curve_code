pro plot_rad_vs_wavl,psplot=psplot,showstarspec=showstarspec,$
                     nbins=nbins,custfile=custfile,$
                     showtheospec=showtheospec,choosefile=choosefile,$
                     totsets=totsets,wavnum=wavnum,custXrange=custXrange,$
                     showOptical=showOptical,custYrange=custYrange,$
                     powerErr=powerErr,multErr=multErr,medianbin=medianbin,$
                     nolit=nolit,showtext=showtext,depthkep=depthkep,$
                     kepMORIS=kepMORIS,phot=phot,custcharS=custcharS,$
                     asymmetric=asymmetric,$
                     rightleg=rightleg,bottomleg=bottomleg,$
                     filterCurveColor=filterCurveColor,$
                     prevChoices=prevChoices
;;psplot -- saves a postscript plot
;;showstarspec -- shows a star spectrum on the same plot
;;nbins -- number of points bo bin in Rp/R*
;;custfile -- chooses a custom radius vs wavelength file
;;showtheospec -- shows a theoretical transmission spectrum
;;choosefile -- choose a file from the radius_vs_wavelength directory
;;totsets -- optional keyword to specify the number of total sets of
;;           data to over-plot (eliminates the old additionalfile keyword)
;;wavnum -- changes the units on wavelength to wavenumber
;;custXrange -- set a custom range for the plot instead of defualt
;;showOptical -- show the optical transit radius
;;custYrange -- set a custom range fot eh plot instead of default
;;powerErr -- takes all radius errors to a specified power
;;multErr -- takes all radius errors and multiples by a specified constant
;;medianbin -- bins the points by the median value instead of the
;;             weighted average
;;nolit -- no literature value for the planet radius
;;showtext -- explanatory text
;;depthkep - labels the Y axis as Kepler transit depth instead of Rp/R*
;;kepMORIS - shows the approx MORIS transit depth
;; phot - shows the zprime data as different
;; custcharS -- set a custom character size (Terry thought my numbers
;;              were too big)
;; asymmetric -- show asymmetric error bars
;; rightleg/bottomleg -- move the legend to the right/bottom
;; filterCurveColor -- lets you specify the filter curve color
;; prevChoices - Use the previous choices for rad_vs_wavl files and
;;               legend labels

  if keyword_set(showstar) then !x.margin = [9,9] else !x.margin=[9,3]
  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/rad_vs_wavl'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=10, ysize=7,decomposed=1,/color
  endif

  ;; restore the star spectrum to show where telluric & star features
  ;; might be
  restore,'data/specdata.sav'

  if keyword_set(prevChoices) then restore,'param_input/rad_file_choices.sav'

  ;; read in the radius versus wavelength file
  case 1 of
     keyword_set(custfile): radfile=custfile
     keyword_set(prevChoices): radfile = prevRadFArr[0]
     keyword_set(choosefile): begin
        radfile = choose_file(searchDir='radius_vs_wavelength',$
                              filetype='.txt')
        ;; Search the radius_vs_wavlength directory for .txt files
     end
     else: radfile='radius_vs_wavelength/radius_vs_wavl.txt'
  endcase


  if keyword_set(asymmetric) then begin
     readcol,radfile,wavl,wavlsize,rad,rade,radep,radem,skipline=1,format='(F,F,F,F)'
  endif else begin
     readcol,radfile,wavl,wavlsize,rad,rade,skipline=1,format='(F,F,F)'
     ;; Make the +/- bars the same
     radep = rade
     radem = rade
  endelse

  if n_elements(powerErr) NE 0 then rade = rade^(powerErr)
  if n_elements(multErr) NE 0 then rade = rade * multErr

  if keyword_set(showstarspec) then ytempstyle=8 else ytempstyle=1

  binsizes = wavlsize

  if n_elements(nbins) NE 0 then begin
     ;; if asked to, bin the Rp/R*
     nnewwavl = ceil(float(n_elements(wavl))/float(nbins))
     newwavl = fltarr(nnewwavl)
     newrad = fltarr(nnewwavl)
     newrade = fltarr(nnewwavl)
     newbinsizes = fltarr(nnewwavl)
     for i=0l,nnewwavl-1l do begin
        subInd = lindgen(nbins)+i*nbins
        ;; check that you haven't gone over the number of points
        allowedsubpts = where(subInd LE n_elements(wavl) -1l)
        allowedpts = subInd[allowedsubpts]
        newwavl[i] = mean(wavl[allowedpts])
        weights = 1E/rade[allowedpts]^2
        if keyword_set(medianbin) then begin
           newrad[i] = median(rad[allowedpts])
        endif else begin
           newrad[i] = total(weights *rad[allowedpts])/total(weights)
        endelse
        newbinsizes[i] = binsizes[i] * float(n_elements(allowedpts))
        newrade[i] = sqrt(1E /total(weights))
     endfor
     wavl = newwavl
     rad = newrad
     rade = newrade
     binsizes = newbinsizes
     ;; Save the binned wavelength file
     forprint,wavl,binsizes,rad,rade,$
              textout='radius_vs_wavelength/binned_rp_rs.txt',$
              comment='# Wavelength(um)  Bin size   Rp/R*    Rp/R* Error'
  endif else begin
  endelse

  myplotsym = 0 ; circle
  myplotfill = 1 ; fill
  myplotsymSize = 0.5 ;; half size
  plotsym,myplotsym,myplotsymSize,fill=myplotfill 

  wavlwidth = binsizes/2E
  if keyword_set(wavnum) then begin
     myxtitle = 'Wave Number (cm!E-1!N)'
     wavlwidth = 1E4 / wavl^2 * wavlwidth
     wavl = 1E4 / (wavl)
     custxrange = [4000,11500]
  endif else begin
     myxtitle='Wavelength ('+cgGreek('mu')+'m)'
  endelse

  if n_elements(depthkep) NE 0 then begin
     myYtitle='Transit Depth / Mean Kepler'
     mylinestyle=1
     wavlwidth = wavlwidth * 0E
     ;; change the default ranges
     if n_elements(custYrange) EQ 0 then custYrange = [-3,9]
     if n_elements(custXrange) EQ 0 then custxrange=[0.5,2.5]
  endif else begin
     myYtitle = 'R!Dp!N/R!D*!N'
     mylinestyle=0
     if n_elements(custXrange) EQ 0 then custxrange=[0.8,2.55]
     if n_elements(custYrange) EQ 0 then custYrange=[0.12,0.17]
  endelse

  if n_elements(custcharS) EQ 0 then custcharS = 1E

  plot,wavl,rad,$
       xtitle=myxtitle,$
       ytitle=myYtitle,$
       ystyle=ytempstyle,xstyle=1,xrange=custxrange,$
       yrange=custYrange,/nodata,charsize=custcharS

  if keyword_set(phot) then begin
     oploterror,[wavl[0]],[rad[0]],[0],[rade[0]],thick=4,$
                hatlength=!D.X_VSIZE / 30,errstyle=0,psym=8,symsize=1.2
     nwavs = n_elements(wavl)
     oploterror,wavl[1:nwavs-1l],rad[1:nwavs-1l],wavlwidth[1:nwavs-1],rade[1:nwavs-1],$
                psym=8,thick=2,linestyle=mylinestyle
  endif else begin
     oploterror,wavl,rad,wavlwidth,radep,psym=8,thick=2,/hibar
     oploterror,wavl,rad,wavlwidth,radem,psym=8,thick=2,/lobar
  endelse
  if keyword_set(depthkep) then oplot,wavl,rad,thick=2,linestyle=0

;                color=mycol('yellow') 
  
  ;; As in Gibson et al. 2012, show 3 scale heights around the
  ;; adopted Rp/R* from Jacob Bean et al. 2012
  scaleH = 0.00115E
;  scaleH = 0.00115E * 2E
  if not keyword_set(nolit) and not keyword_set(depthkep) then begin
     plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433],color=mycol(['red'])
     plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]+3E*scaleH,color=mycol(['red']),linestyle=2
     plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]-3E*scaleH,color=mycol(['red']),linestyle=2
  endif
        
     
  prevXrange = !x.crange
  if keyword_set(showtheospec) then begin ;; show the theoretical transmission spectrum
;     readcol,'../models/fortney_g10mps_2500K_isothermal.csv',theowav,theorad,$
;             skipline=6,format='(F,F)'
;      readcol,'../models/transit_models/transit_t2500g10_noTiO.dat',theowav,theorad,$
;              format='(F,F)'
;      readcol,'../models/transit_models/lambda_R_P_iso_g10_2500.dat',theowav,theorad,$
;              format='(F,F)'
;     mult2 = 0.10418 ;; found from the minimum chi-squared
;     mult2 = 1.54369E-6 ;; found from the minimum chi-squared
;     mult2 = 1.45734E-6 ;; found from the minimum chi-squared

;     theorad = smooth(theorad,5)
     restore,'data/binned_final.sav'
     modcolor = mycol(['blue','dgreen'])
     positionMultX = [1.0,0.63] ;; position multipliers
     positionMultY = [1.018,0.99] ;; position multipliers

     for i=0l,nmod-1l do begin
        ;; Full resolution
        xsmooth = fullRes.(i * 2l)
        ySmooth = specsmooth(xsmooth,fullRes.(i * 2l+1l),100)
        oplot,xsmooth,ysmooth,color=modcolor[i]
        ;; Binned
        plotsym,myplotsym,myplotsymSize,fill=0
        oplot,binnedWav,binnedValues[*,i],psym=8,color=modcolor[i],symsize=1
        xyouts,binnedWav[0] * positionMultX[i],binnedValues[0,i] * positionMultY[i],$
               modName[i],charsize=0.7,$
               color=modColor[i]
     endfor

;     ntheo=n_elements(theorad)
     ;; Bin model over wavelenght ranges
;     binModel = avg_series(theowav,theorad *mult2,fltarr(ntheo)+0.2E,wavl-wavlwidth,wavlwidth * 2E,weighted=0)

     if keyword_set(showOptical) then begin
        ;; Show the Bean 2009 result if asked to
        CorRad = 0.1433 ;;Rp/R*
        CorErr = 0.0010
        CorWav = 0.65 ;; microns, approximately
        CorWid = 0.20 ;; microns, approx
;        binModel2 = avg_series(theowav,theorad * mult2,fltarr(ntheo)+0.2E,CorWav-CorWid/2E,CorWid,weighted=0)
;        oplot,[CorWav],[binModel2],psym=2,color=mycol('blue'),symsize=2
        if keyword_set(phot) then begin
           CorWid = 0E ;; instead, we'll show the filter curve
           myCorSymbol = 4
        endif else myCorSymbol = 3
        oploterror,CorWav,CorRad,CorWid,CorErr,color=mycol('red'),thick=2,psym=myCorSymbol,symsize=0.7
     endif
  endif

  if keyword_set(kepMORIS) then begin
     morisRad = 0.5
     morisRadE = 0
     morisWav = 0.626
     morisWid = 0.1342
     oploterror,morisWav,morisRad,morisWid,morisRadE,color=mycol('red'),psym=3,thick=2
  endif

  ;; if undefined, only show one set of data
  if n_elements(totsets) EQ 0 then totsets=1
  radFarrSave = replicate('',totsets)
  radFarrSave[0] = radfile


  if totsets GT 1l then begin
     legnamearr = strarr(totsets)
     print,'Name for file '+radfile+' ?'
     tempnm=''
     if keyword_set(prevChoices) then begin
        tempnm = prevLegNmArr[0]
     endif else read,tempnm,format='(A)'
     legnamearr[0l] = tempnm
     colorlist = [!p.color,mycol(['orange','purple','blue'])]
     ncolchoices = n_elements(colorlist)
     colorchoices = colorlist[lindgen(totsets) mod ncolchoices]
  endif
     
  for i=2l,totsets do begin
     ;; If asked to, overplot another Rad/vs wavlength file
     print,'Choose Additional file ',strtrim(i-1l,1)
     if keyword_set(prevChoices) then begin
        file2 = prevRadFarr[i-1]
     endif else file2 = choose_file(searchDir='radius_vs_wavelength',filetype='.txt')
     readcol,file2,wavl2,wavl2size,rad2,rade2,skipline=1,format='(F,F,F)'
     ;; find the bin width
     wavlwidth2 = wavl2size/2E
     if keyword_set(wavnum) then begin
        wavlwidth = 1E4 / wavl^2 * wavlwidth
        wavl = 1E4 / (wavl)
     endif

     if keyword_set(depthkep) then wavlwidth2 = wavlwidth2 * 0E

     if keyword_set(phot) then begin
        oploterror,[wavl2[0]],[rad2[0]],[0],[rade2[0]],thick=4,psym=8,$
                   color=colorchoices[i-1l],errstyle=0,hatlength=!D.X_VSIZE / 30,$
                   symsize=1.2
        nwavs = n_elements(wavl2)
        oploterror,wavl2[1:nwavs-1l],rad2[1:nwavs-1l],wavlwidth2[1:nwavs-1],rade2[1:nwavs-1],$
                   psym=8,thick=2,color=colorchoices[i-1l]
     endif else begin
        oploterror,wavl2,rad2,wavlwidth2,rade2,psym=8,thick=2,color=colorchoices[i-1l]
     endelse

     if keyword_set(depthkep) then oplot,wavl2,rad2,thick=2,linestyle=0,$
                                         color=colorchoices[i-1l]

     print,'Legend Name for data from file '+file2+' ?'
     tempnm = ''
     if keyword_set(prevChoices) then begin
        tempnm = prevLegNmArr[i-1]
     endif else read,tempnm,format='(A)'
     legnamearr[i-1l] = tempnm
     radFarrSave[i-1] = file2
  endfor

  ;; Save all the file name choices and Legend name choices in case
  ;; you want to use them again
  if totsets GT 1 then prevLegNmArr = legnamearr else prevLegNmArr = ''
  prevRadFarr = radFarrSave
  save,prevLegNmArr,prevRadFarr,filename='param_input/rad_file_choices.sav'

  if keyword_set(psplot) then legcharsize = 0.75 else legcharsize=1

  if totsets GT 1l then begin
     if n_elements(rightleg) EQ 0 then begin
        if keyword_set(depthk) then rightleg=1 else rightleg=0
     endif
     if n_elements(bottomleg) EQ 0 then bottomleg = 0
     al_legend,legnamearr,psym=8+lonarr(totsets),color=colorchoices,/clear,$
               right=rightleg,bottom=bottomleg,charsize=legcharsize,symsize=1.05 + lonarr(totsets)

  endif

  if keyword_set(showtext) then begin
     plotsym,myplotsym,myplotsymSize,fill=1
     al_legend,['Bean 2009','This Work','Binned Model Value'],psym=[4,8,8],/bottom,charsize=legcharsize,$
            color=[mycol('white'),!P.color,mycol('white')],/right,fill=1,symsize=[0.7,1.0,1.0],$
               thick=[2,1,1],textcolor=mycol(['white','white','white'])
     plotsym,myplotsym,myplotsymSize,fill=0
     al_legend,['Bean 2009','This Work','Binned Model Value'],psym=[4,8,8],/bottom,charsize=legcharsize,$
            color=[mycol('red'),!P.color,!P.color],/right,fill=1,symsize=[0.7,1.0,1.0],$
               thick=[2,1,1]
  endif

  if keyword_set(phot) then begin
     readcol,'../calculations/zprime_transmission/zprime_response.txt.csv',skipline=1,$
             wavel,trans
     if n_elements(filterCurveColor) EQ 0 then filterCurveColor='dgreen'
     oplot,wavel,trans / max(trans) * 0.1E * (!y.crange[1] - !y.crange[0]) + !y.crange[0],$
           color=mycol(filterCurveColor)
     if keyword_set(showOptical) then begin
        readcol,'../corot_data/filter_curve/filter_curve_CoRoT.txt',skipline=1,$
                CoRoTtransWav,CoRoTtrans
        CoRoTtransWav = CoRoTtransWav * 1E-3 ;; convert from nm to microns
        oplot,CoRoTtransWav,CoRoTtrans / max(CoRoTtrans) * 0.1E * (!y.crange[1] - !y.crange[0]) + !y.crange[0],$
              color=mycol(filterCurveColor),linestyle=2
        
     endif
  endif

  if keyword_set(showstarspec) then begin
     ;; plot the source spectrum
     plot,lamgrid,flgrid(*,0,1),/noerase,xrange=prevXrange,ystyle=5,xstyle=1,$
          /nodata;yrange=[-6E5,6E5],/nodata
     oplot,lamgrid,flgrid(*,0,1),color=mycol('blue')
     axis,yaxis=1,yrange=!y.crange,color=mycol('blue'),/ystyle,$
          ytitle='Raw Source Flux (DN)'
     
     !x.margin = [10.0,3.0]
  endif

  if keyword_set(depthkep) then begin
     oplot,!x.crange,[1,1],linestyle=2,color=mycol('red')
  endif

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 300% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
  !x.margin = [15.0,4.0]

end  
  
