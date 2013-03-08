pro plot_sim,range=range,photcompar=photcompar,$
            gphot=gphot,$
            centcomparX=centcomparX,$
            centcomparY=centcomparY,$
            fwhmcompar=fwhmcompar,$
            xfwhmcompar=xfwhmcompar,$
            yfwhmcompar=yfwhmcompar,$
            thetacompar=thetacompar,$
            reftype=reftype,$
            testmod=testmod,$
            weighted=weighted,$
            simplegreen=simplegreen,$
            nosigreject=nosigreject,$
            comptype=comptype,$
             fullspecrange=fullspecrange,$
             psplot=psplot
;; This script is designed to plot spectral data and guider data
;; simultaneously

;;range -- allows you to specify the plot range
;;photcompar -- will overplot photometry info
;;centcomparX -- will overplot the y centroids
;;centcomparY -- will overplot the y centroids
;;fwhmcompar -- will overplot the average FWHM
;;Xfwhmcompar -- will overplot the x FWHM
;;Yfwhmcompar -- will overplot the y FWHM
;;thetacompar -- will overlpot the angle with respect to the horizontal
;;reftype -- sets whether you're talking about the source star or the reference stars
;;testmod -- will test out a  model between spectral shift and photometry data
;;weighted -- (default, true) uses the weighted average of photometry
;;            data over spectral exposures
;;trycorrect -- attempts to divide the shifts by a previousy
;;              calculated model shift
;;orderl -- sets the order to look at
;;simplegreen -- makes the plots less busy by making all photometry green
;;nosigreject -- turns of the rejection algorithm when binning
;;               photometric data
;;comptype -- chooses the type of spectral data to compare to (default
;;            is wavelength shift)
;;fullspecrange -- show all the spectral data in the correlation plot
;;                 - not just the middle of it
;;psplot -- saves postscript plots instead of making X windows plots

if keyword_set(psplot) then begin
   set_plot,'ps'
   !p.font=0
endif else begin
   window,1,ypos=20
endelse

UTxrange = [10000,20000]
;UTPosB = [8.606,8.622] ;UT times where the nod position is B

;common phD,t,explist
; Get the photometry data
;; Get the filename parameters
readcol,'param_input/geomap_filenm.txt',filedescrip,filen,$
        format='(A,A)'

;; Read in the geomap data
readcol,filen[0],/silent,$
        filenl,xshift,xrms,yshift,yrms,xrot,yrot,$
        format='(A,F,F,F,F,F,F)',skipline=1

;; get the UT times, JD times & exptimes of the photometry
readcol,filen[1],/silent,$
        filenl2,gexpt,gtime,date,JD,$
        format='(A,F,A,A,D)'

;; Put the photometry info into a structure
photLength = n_elements(yshift)
stnumarr = intarr(photLength) + 1
phD = create_struct('stnum',stnumarr,$ ;; one star since the geomap output is for 1 star
                    'UTSHUT',JD,$
                    'y',yshift,$
                    'x',xshift,$
                    'theta',xrot)



;readcol,'bx_wavelengths.txt',bx,wavl,skipline=1,$
;  format='(I,F)'
;nbox = n_elements(bx)
;;readcol,'fl_time_ser.txt',t,flA,flB,flC,flD,FlE,format='(D,F,F,F,F,F,F,F,F,F)'
;nrows = file_lines('fl_time_ser.txt') ; number of rows
;nrows -= 1l ;; subtract out for the header
;;nbox = file_lines('ds9_boxes.reg') ; number of boxes extracted
;header = strarr(1)

if keyword_set(simplegreen) then begin
   colorOpt=mycol(['green'])
   psymOpt = [2]
endif else begin
   colorOpt = mycol(['white','red','green','lblue','magenta','yellow',$
                     'pink','orange','blue'])
   psymOpt = [1,2,4,5,6,7] ;; plot symbols that are actually readable
endelse

ncol = n_elements(colorOpt)
nsym = n_elements(psymOpt)


lineOpt = [0,2,3,4,5]


readcol,'../IRTF_UT2012Jan04/raw/spectrograph_ut_and_exp_list.txt',$
        specfilen,t,explist,format='(A,D,F)'
nrows = n_elements(t) ;; number of spectra files

assert,nrows,'=',n_elements(explist),$
       'Exposure time list mismatches list of UTSHUT times'
;
;if n_elements(range) EQ 0 then yRan=[0.65,1.2] else begin
;   yRan = range
;endelse
;
;window,1,xsize=700,ysize=500,ypos=20
;
!x.margin = [13,14]
;

;; Read in the spectral data
restore,'data/specdata.sav'
bingridmiddle = bingrid + binsizes/2E

wavindex = 4
wavname = string(bingridmiddle[wavindex],format='(F4.2)')
if n_elements(comptype) GT 0 then begin
   case comptype of
      'ratio': begin
         wsha = transpose(binfl[wavindex,*])
         compname='Flux Ratio at '+wavname+' um'
         compunits=''
         fcompname=strtrim(wavname,1)+'_ratio'
         fpref = strtrim(wavname,1)+'_ratio'
      end
      'host': begin
         wsha = transpose(binind[wavindex,0,*])
         compname='Host Star Flux at '+wavname+' um'
         compunits='(DN)'
         fcompname='host_'+strtrim(wavname,1)+'_fl'
         fpref = fcompname
      end
      'refstar':begin
         wsha = transpose(binind[wavindex,1,*])
         compname='Ref Star Flux at '+wavname+' um'
         compunits='(DN)'
         fcompname='ref_'+strtrim(wavname,1)+'_fl'
         fpref = fcompname
      end
      else: begin
         print,'Unrecognized Comparison Type'
         return
      end
   endcase
endif else begin
   wsha = transpose(binfl[wavindex,*])
   compname='Flux Ratio at '+wavname+' um'
   compunits=''
   fcompname='mean_'+strtrim(wavname,1)+'_fl'
   fpref = 'mean_'+strtrim(wavname,1)+'_fl'
endelse



case 1 of 
   keyword_set(centcomparY): fsuffix = 'y'
   keyword_set(centcomparX): fsuffix = 'x'
   keyword_set(fwhmcompar): fsuffix = 'fwhm_av'
   keyword_set(xfwhmcompar): fsuffix = 'fwhm_x'
   keyword_set(yfwhmcompar): fsuffix = 'fwhm_y'
   keyword_set(gphot): fsuffix = 'gphot'
   keyword_set(photcompar): fsuffix = 'apphot'
   keyword_set(testmod): fsuffix = 'mod'
   keyword_set(thetacompar): fsuffix = 'theta'
   else: fsuffix = ''
endcase


if keyword_set(psplot) then begin
   plotfn1 = 'plots/specphot_cors/time_series_'+fsuffix
   device,encapsulated=1, /helvetica,$
          filename=plotfn1+'.eps'
   device,xsize=14, ysize=10,decomposed=1,/color
endif

;; Turn all times into hours after an arbitrary time
JDstart = date_conv('2012-01-04T08:23:00','JULIAN')
t = (t - JDstart)*24D*3600D
phD.UTSHUT = (phD.UTSHUT - JDstart)*24D*3600D


assert,nrows,'=',n_elements(wsha),$
       'Number of files mismatch between shift and exposure times'

plot,t,wsha,ytitle=compname+' '+compunits,xtitle='t (sec)',$
     ystyle=8,xrange=UTxrange,/xstyle,/nodata

;; plot the wavelength shift as a bar to show the exposure time
for j=0l,nrows-1l do begin
   plots,[t[j],t[j]+explist[j]],[wsha[j],wsha[j]]
endfor


;exptime= 60D / (3600D) ; exposure time is 60 seconds for HD8991
;stdarr = dblarr(nbox)
;errarr = dblarr(nbox)
;for j=0l,nbox-1l do begin
;;   oplot,t,nfl[j,*],color=colorarr[j],linestyle=linarr[j]
;    for k=0l,8l - 1l do begin
;        plots,[t[k],t[k]+exptime],[nfl[j,k],nfl[j,k]],$
;          color=colorarr[j],linestyle=linarr[j]
;    endfor
;
;   stdarr[j]=stdev(nfl[j,*])
;   errarr[j]=sqrt(median(fl[j,*]))/median(fl[j,*])
;   print,'Box ',strtrim(j,1),' stdev/counting noise = ',stdarr[j]/errarr[j]
;endfor
;
;legend,'Box '+strtrim(lindgen(nbox),1),$
;       color=colorarr,$
;       linestyle=linarr,/right
;
;forprint,lindgen(nbox),stdarr,errarr,textout='standard_dev_array.txt',comment='Box #   Standard Deviation of normalized flux   Stdev/counting noise'
;



if keyword_set(centcomparY) OR $
   keyword_set(centcomparX) OR $
   keyword_set(fwhmcompar) OR $
   keyword_set(xfwhmcompar) OR $
   keyword_set(yfwhmcompar) OR $
   keyword_set(gphot) OR $
   keyword_set(photcompar) OR $
   keyword_set(testmod) OR $ 
   keyword_set(thetacompar) $
then begin

;restore,"flt_ph.sav"
nsrc = max(phD.stnum)

;if n_elements(reftype) NE 0 then begin ;; choose whether to compare to reference stars or source
;    case reftype of
;        'source': refrange = [1,1]
;        'references': refrange = [2,nsrc]
;        else: refrange = [2,nsrc]
;    endcase
;    print,"reftype = ",reftype
;endif else refrange = [2,nsrc]
refrange = [1,nsrc]

;; start an array for weighted averages
ysub = dblarr(n_elements(phD.UTSHUT))*!values.d_nan

;; make an array for the colors
colorarr = colorOpt[lindgen(nsrc+1l) mod ncol]
 psymarr =  psymOpt[lindgen(nsrc+1l) mod nsym]

if keyword_set(photcompar) then PHOTarr = phD.phot else begin
    if keyword_set(gphot) then PHOTarr = phD.gphot
endelse
if n_elements(PHOTarr) NE 0 then begin
;    plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
   stind = where(phD.stnum EQ refrange[0])
   avg = mean(PHOTarr[stind])
   topY = 1.2*max(PHOTarr[stind],/nan)/avg
   botY = 0.5*min(PHOTarr[stind],/nan)/avg
   plot,phD.UTSHUT,PHOTarr[stind]/avg,$
        /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
        yrange=[botY,topY]
;      yrange=[0.1,2.5]
    ynam = 'Normalized Phot'
    fsuffix = 'phot'

    for i=refrange[0],refrange[1] do begin
        ;;Normalize and plot
        stind = where(phD.stnum EQ i)
;        if n_elements(stind) GT 30 then begin
        if stind NE [-1] then begin
            avg = mean(PHOTarr[stind])
            oplot,phD.UTSHUT[stind],PHOTarr[stind]/avg,color=colorarr[i],$
                  psym=psymarr[i]
            ysub[stind] = PHOTarr[stind]/avg
        endif
    endfor

endif

if keyword_set(centcomparY) then begin
    plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
      /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
      yrange=[-5,5]
;    AXIS, YAXIS=1, YRANGE=!y.crange, ytitle='Y centroid shift
;    (px)',color=green,/ystyle
    ynam='Y centroid shift (px)'
    fsuffix = 'centY'
    for i=refrange[0],refrange[1] do begin
        ;;Normalize and plot
        stind = where(phD.stnum EQ i)
        if stind NE [-1] then begin
            avg = mean(phD.y[stind])
            dy = phD.y[stind] - avg
            oplot,phD.UTSHUT[stind],dy,color=colorarr[i],$
                  psym=psymarr[i]
            ysub[stind] = (phD.y[stind]-avg)
        endif
    endfor
endif

if keyword_set(centcomparX) then begin
    plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
      /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
         yrange=[-10,10]

;         yrange=[-5,5]
;         yrange=[-4,3.5]
;    AXIS, YAXIS=1, YRANGE=!y.crange, ytitle='Y centroid shift
;    (px)',color=green,/ystyle
    ynam='X centroid shift (px)'
    fsuffix = 'centX'
    for i=refrange[0],refrange[1] do begin
        ;;Normalize and plot
        stind = where(phD.stnum EQ i)
;        bpos = where(phD.UTSHUT[stind] GT UTPosB[0] and phD.UTSHUT[stind] LE UTPosB[1])
        if stind NE [-1] then begin
           ;; find average at the B position
           avg = mean(phD.x[stind])
           oplot,phD.UTSHUT[stind],(phD.x[stind]-avg),color=colorarr[i],$
                 psym=psymarr[i]
           ysub[stind] = (phD.x[stind]-avg)
        endif
    endfor
endif

if keyword_set(fwhmcompar) then begin
    plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
      /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
      yrange=[0,7]
    ynam = 'Average FWHM (px)'
    fsuffix='av_fwhm'

    ;; find average of major and minor FWHM
    mmFWHM = (phD.x_FWHM + phD.y_FWHM)/2E

    ;; multiply by the average FWHM for all stars
    fullFWHMavg = mean(mmFWHM,/nan)

    for i=refrange[0],refrange[1] do begin
        ;;Normalize and plot

        stind = where(phD.stnum EQ i)
        if stind NE [-1] then begin
            avg = mean(mmFWHM[stind])
            ysub[stind] = mmFWHM[stind] / avg * fullFWHMavg
;;            ysub[stind] = mmFWHM[stind]
            oplot,phD.UTSHUT[stind],ysub[stind],color=colorarr[i],$
                  psym=psymarr[i]
        endif
    endfor

endif

if keyword_set(Xfwhmcompar) then begin
    plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
      /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
      yrange=[0,7]
    ynam = 'X FWHM (px)'
    fsuffix='Xfwhm'
    ;; multiply by the average FWHM for all stars
    fullFWHMavg = mean(phD.x_FWHM,/nan)

    for i=refrange[0],refrange[1] do begin
        ;;Normalize and plot

        stind = where(phD.stnum EQ i)
        if stind NE [-1] then begin
            avg = mean(phD.x_FWHM[stind])
            ysub[stind] = phD.x_FWHM[stind] / avg * fullFWHMavg
;;            ysub[stind] = phD.x_FWHM[stind]
            oplot,phD.UTSHUT[stind],ysub[stind],color=colorarr[i],$
                  psym=psymarr[i]
        endif
    endfor

endif

if keyword_set(Yfwhmcompar) then begin
    plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
      /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
      yrange=[0,7]
    ynam = 'Y FWHM (px)'
    fsuffix='Y_fwhm'
    ;; multiply by the average FWHM for all stars
    fullFWHMavg = mean(phD.y_FWHM,/nan)

    for i=refrange[0],refrange[1] do begin
        ;;Normalize and plot

        stind = where(phD.stnum EQ i)
        if stind NE [-1] then begin
            avg = mean(phD.y_FWHM[stind])
;;            ysub[stind] = phD.x_FWHM[stind] / avg * fullFWHMavg
            ysub[stind] = phD.y_FWHM[stind]
            oplot,phD.UTSHUT[stind],ysub[stind],color=colorarr[i],$
                  psym=psymarr[i]
        endif
    endfor

endif

if keyword_set(thetacompar) then begin

   ;;make the rotated theta aray
   mytheta = fltarr(n_elements(phD.theta))
   lowp = where(phD.theta LT 20E,complement=highp)
   mytheta = fltarr(n_elements(phD.theta))
   mytheta[lowp] = 180E + phD.theta[lowp]
   if highp NE [-1] then begin
      mytheta[highp] = phD.theta[highp]
   endif

   plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
      /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
      yrange=[-50,300]
    ynam = 'Adjusted theta (deg)'
    fsuffix='theta'
    ;; multiply by the average FWHM for all stars
;    fullFWHMavg = mean(phD.x_FWHM,/nan)

    for i=refrange[0],refrange[1] do begin
        ;;Normalize and plot

        stind = where(phD.stnum EQ i)
        if stind NE [-1] then begin
;;            avg = mean(phD.x_FWHM[stind])
;;            ysub[stind] = phD.x_FWHM[stind] / avg * fullFWHMavg
            ysub[stind] = mytheta[stind]
            oplot,phD.UTSHUT[stind],ysub[stind],color=colorarr[i],$
                  psym=psymarr[i]
        endif
    endfor

endif

if keyword_set(testmod) then begin
   if comptype EQ 'slope' then begin
      plotmodrange=[-3E3,3E3]
   endif else plotmodrange = [-1.5E-4,1.5E-4]
    plot,phD.UTSHUT,dblarr(n_elements(phD.UTSHUT))+1D,$
      /nodata,/noerase,xrange=UTxrange,/xstyle,ystyle=5,$
      yrange=plotmodrange

    ynam = 'Model '+compname
    fsuffix='modeltest'

;    nPhot = dblarr(n_elements(phD.UTSHUT))*!values.d_nan
;    nSourcePhot = dblarr(n_elements(phD.UTSHUT))*!values.d_nan
    FWHM = dblarr(n_elements(phD.UTSHUT))*!values.d_nan
;    FWHMY = dblarr(n_elements(phD.UTSHUT))*!values.d_nan
    Ypos = dblarr(n_elements(phD.UTSHUT))*!values.d_nan
    Xpos = dblarr(n_elements(phD.UTSHUT))*!values.d_nan

    ;;Find the average FWHM for all stars
    fullavFWHM = mean((phD.x_FWHM + phD.y_FWHM)/2E,/nan)

    ;; initialize arrays for the photometry
    for i=2,nsrc do begin
       stind = where(phD.stnum EQ i) ;; find the indices of the star
       bpos = where(phD.UTSHUT[stind] GT UTPosB[0] and phD.UTSHUT[stind] LE UTPosB[1])
       if stind NE [-1] and bpos NE [-1] then begin
          ;;nPhot[stind] = phD.phot[stind] / mean(phD.phot[stind]) ;;normalized photometry of reference stars
;          FWHMX[stind] = 
          FWHM[stind] = (phD.x_FWHM[stind] + phD.y_FWHM[stind])/2E
          FWHM[stind] = FWHM[stind] / mean(FWHM[stind]) * fullavFWHM;; reNormalize

          Ypos[stind] = phD.y[stind] - mean(phD.y[stind])
          avgX = mean(phD.x[stind[bpos]])
          Xpos[stind] = phD.x[stind] - avgX
          
;          ysub[stind] = -2E-4 *slit_trans_mom(Ypos[stind]-0.2E,FWHMY[stind],4.27E/1.3E)
;          oplot,phD.UTSHUT[stind],ysub[stind],psym=2,color=mycol('green')

       endif
    endfor


    ;; find the best fit model assuming a moment of a gaussian
    ;; convolved with a rectangular slit

    nphot = n_elements(phD.UTSHUT) ;; number of photometric extractions
    inarr = fltarr(nphot,8);; input array for the model
    inarr[*,0] = Ypos
    inarr[*,1] = FWHM
    inarr[*,2] = phD.UTSHUT
    inarr[*,3] = phD.SNR_PH
    inarr[*,4] = fltarr(nphot) + float(nrows) ;; tells how many spectra there are
    inarr[0l:nrows-1l,5] = t
    inarr[nrows:nphot-1l,5] = 0E
    inarr[0l:nrows-1l,6] = explist
    inarr[nrows:nphot-1l,6] = 0E
    inarr[*,7] = Xpos

    ;; parameter info array
    pi = replicate({fixed:0,limited:[0,0],limits:[0E,0E]},3)
    case testmod of
       'GausM':begin ;; moment of slit transmission model
          expr = 'model_moment_eval(X,P[0],P[1],P[2])'
          parnames = ['amplitude','y offset','y slit width']
          start = [-2E-4,0E,4.27E]
       end
       'GausS':begin ;; model of slit transmission models' slope
          expr = 'model_slope_eval(X,P[0],P[1],P[2])'
          parnames = ['amplitude','y offset','y slit width']
          start = [1E6,0.1E,4.27E]
          pi(2).limited(0) = 1
          pi(2).limits(0) = 0.5E ;; assert that the slit width must be greater than 2 pixels
          pi(2).limited(1) = 1
          pi(2).limits(1) = 9E ;; assert that the slit width must be less than 8 pixels
       end
       'GausT':begin ;; model of slit transmission
          expr = 'model_trans_eval(X,P[0],P[1],P[2])'
          parnames = ['amplitude','y offset','y slit width']
          start = [1E6,0.1E,4.27E]
          pi(2).limited(0) = 1
          pi(2).limits(0) = 0.5E ;; assert that the slit width must be greater than 2 pixels
          pi(2).limited(1) = 1
          pi(2).limits(1) = 9E ;; assert that the slit width must be less than 8 pixels
       end
       else: begin
          expr = 'lin_wav_mod(X,P[0],P[1],P[2])'
          parnames = ['A','y0','B']
          start = [-2E-4,0E,4.27E]
          testmod='linMod'
       endelse
    endcase

    
    ;; stick in a uniform array for errors
    wshaerr = mean(wsha)/20E
    if comptype EQ 'slope' then wshaerr += 0.001E ;; make sure errors aren't zero
    result=mpfitexpr(expr,inarr,wsha,wshaerr,start,PARINFO=pi,perr=punct)

    openw,1,'shift_info/model_params_'+fpref+'.txt'
    printf,1,'Model for '
    for i=0,n_elements(parnames)-1 do begin
       printf,1,parnames[i],' = ',result[i]
    endfor
    printf,1,'See model_moment_eval.pro for model'
    close,1

    case testmod of
;       'GausM':ysub = result[0] *slit_trans_mom(Ypos-result[1],FWHM,result[2])
;       'GausS':ysub = model_slope_eval(inarr,result[0],result[1],result[2],/nobin)
;       'GausT':ysub = model_trans_eval(inarr,result[0],result[1],result[2],/nobin)
       else: ysub = result[0] * (Ypos - result[1]) + result[2] * Xpos
    endcase

    oplot,phD.UTSHUT,ysub,psym=2,color=mycol('blue')
;    stind = where(phD.stnum EQ 1) ;; indices of source star
;    nSourcePhot[stind] = phD.phot[stind] / mean(phD.phot[stind])

    avgY = avg_series(phD.UTSHUT,ysub,phD.SNR_PH,t,explist,weighted=weighted)
    
    ;; save the predictions
    forprint,wsha,avgY,comment='# Actual shift   Predicted shift (um) '+pref,$
             textout='shift_info/predicted_'+fpref+'s.txt'

 endif




;; show the other variable's axis on the right
AXIS, YAXIS=1, YRANGE=!y.crange,color=mycol('blue'),/ystyle,$
      ytitle=ynam

;; find the weighted average over exposures and plot it

;; Use a uniform Signal to noise raio
SNR=50E + fltarr(photLength)
if keyword_set(nosigreject) then begin
   avgY = avg_series(phD.UTSHUT,ysub,SNR,t,explist,weighted=weighted,eArr=eArr)
endif else begin
   avgY = avg_series(phD.UTSHUT,ysub,SNR,t,explist,weighted=weighted,oreject = 3,eArr=eArr)
endelse




;    avgPhot = avg_series(phD.UTSHUT,nPhot,phD.SNR_PH,t,explist)
;    avgSourcePhot = avg_series(phD.UTSHUT,nSourcePhot,phD.SNR_PH,t,explist)
;    avgFWHMX = avg_series(phD.UTSHUT,FWHMX,phD.SNR_PH,t,explist)
;    avgFWHMY = avg_series(phD.UTSHUT,FWHMY,phD.SNR_PH,t,explist)
;    avgYpos = avg_series(phD.UTSHUT,Ypos,phD.SNR_PH,t,explist)
;    avgXpos = avg_series(phD.UTSHUT,Ypos,phD.SNR_PH,t,explist)
;
;    ;; the slit is 1.0 arcseconds, so that's 4.27 pixels
;    avgY = slit_trans_mom(avgYpos,avgFWHMY,4.27E)
;
;    ;; Save all the averaged data to a structure and a text file
;    avgD = create_struct('t_hrs',t,$
;                         'ref_FWHMx',avgFWHMX,$
;                         'ref_FWHMy',avgFWHMY,$
;                         'ref_y',avgYpos,$
;                         'ref_x',avgXpos)
;;                  'ref_phot',avgPhot,$
;;                  'src_phot',avgSourcePhot,$
;    save,filename='avg_photinfo/avg_data.sav',avgD
;    write_csv,'avg_photinfo/avg_data.csv',avgD,header=tag_names(avgD)


avgT = t + explist/2D
oplot,avgT,avgY,color=mycol('blue'),psym=1,symsize=3,thick=3


;; save the result
if keyword_set(psplot) then begin
   device, /close
   cgPS2PDF,plotfn1+'.eps',$
            /delete_ps
   spawn,'convert -density 160% '+plotfn1+'.pdf '+plotfn1+'.png'
   plotfn2 = 'plots/specphot_cors/'+fpref+'_vs_'+fsuffix
   device,encapsulated=1, /helvetica,$
          filename=plotfn2+'.eps'
   device,xsize=14, ysize=10,decomposed=1,/color
endif else begin
   window,2,xpos=80
endelse
;im = tvrd(true =1)
;if keyword_set(simplegreen) then tssuffix = '_greencol' else begin
;   tssuffix = ''
;endelse
;write_png,'plots/time_series_'+fsuffix+tssuffix+'.png',im



if keyword_set(fullSpecRange) then xdynam = [0,0] else begin
   ;; Show middle the data range (excluding extremes)
   sorty = sort(wsha)
   xlength = n_elements(wsha)
   xlowerL = wsha[sorty[ceil(5E/100E*float(xlength))]] * 0.95
   xUpperL = wsha[sorty[ceil(95E/100E*float(xlength))]] * 1.05
   xdynam = [xlowerL,xUpperL]
endelse

plot,wsha,avgY,psym=1,xtitle=compname+' '+compunits,$
     ytitle='Avg '+ynam,/nodata,ystyle=16,$
     xrange=xdynam
                                ;,yrange=[0.6,1.3]

oploterror,wsha,avgY,fltarr(n_elements(wsha)),eArr,psym=1


;; color by nod position
;nodA1 = where(lindgen(n_elements(wsha)) mod 4 EQ 0)
;nodB1 = where(lindgen(n_elements(wsha)) mod 4 EQ 1)
;nodB2 = where(lindgen(n_elements(wsha)) mod 4 EQ 2)
;nodA2 = where(lindgen(n_elements(wsha)) mod 4 EQ 3)

;oplot,wsha[nodA1],avgY[nodA1],psym=2,symsize=1.5,color=mycol('green')
;oplot,wsha[nodB1],avgY[nodB1],psym=2,symsize=1.5,color=mycol('red')
;oplot,wsha[nodB2],avgY[nodB2],psym=2,symsize=1.5,color=mycol('red')
;oplot,wsha[nodA2],avgY[nodA2],psym=2,symsize=1.5,color=mycol('green')

oplot,wsha,avgY,psym=2,symsize=1.5


if keyword_set(fwhmcompar) OR $
  keyword_set(gphot) OR $
  keyword_set(centcomparX) OR $
  keyword_set(centcomparY) OR $
  keyword_set(thetacompar) OR $
  keyword_set(testmod) OR $
  keyword_set(photcompar) then begin
   
   
   ;;perform a fit of what we called avgY to the average spectral data
   ;;(shift, slope etc)

   ;; First remove all the points that aren't finite
   goodp = where(finite(wsha) EQ 1 and finite(avgY) EQ 1)
   fitY = linfit(wsha[goodp],avgY[goodp])
   
   oplot,wsha[goodp],fitY[0] + fitY[1]*wsha[goodp],color=mycol('blue')

endif

print,"Pearson Correlation Coefficient= ",Correlate(wsha[goodp],avgY[goodp])

;im = tvrd(true=1)
;write_png,'plots/'+fpref+'_vs_'+fsuffix+'.png',im

;chisq = total( (wsha - (fitY[0] + fitY[1]*avgY))^2 )
;print,'Chi-squared for linear fit = ',chisq

endif else begin
;im = tvrd(true =1)
;write_png,'plots/time_series.png',im
endelse

if keyword_set(psplot) then begin
   device,/close
   device,decomposed=0
   cgPS2PDF,plotfn2+'.eps',$
            /delete_ps
   spawn,'convert -density 160% '+plotfn2+'.pdf '+plotfn2+'.png'
   set_plot,'x'
   !p.font=-1
endif

print,'Stdev of '+compname+' '+compunits,stddev(wsha)


end
