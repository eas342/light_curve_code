pro plot_tim_ser,fitrad=fitrad,fitpoly=fitpoly,usepoly=usepoly,npoly=npoly,$
                 fullrange=fullrange,smartbin=smartbin,oneprange=oneprange,$
                 offtranserr=offtranserr,freelimb=freelimb,clarlimb=clarlimb,$
                 psplot=psplot,noreject=noreject,differential=differential,$
                 individual=individual,pngcopy=pngcopy,freeall=freall,fixall=fixall,$
                 timebin=timebin
;; plots the binned data as a time series and can also fit the Rp/R* changes
;; apPlot -- this optional keyword allows one to choose the aperture
;;           to plot
;; fitrad -- this fitting procedure ueses
;;            planet transit information to fit Rp/R* as a
;;            function of wavelength
;; fitpoly -- this fits a polynomial to the data to take out long term
;;           trends
;; usepoly -- this uses the polynomial fit from other sources
;;            (calibrators) on a given source
;; npoly -- allows one to set the number of polynomial parameters used
;; fullrange -- allows for plotting of the full y range
;; smartbin -- use the smart-binned data instead of the binned data
;; oneprange -- sets all the y ranges to be 1%
;; offtranserr -- makes the error equal to the stddev of off-transit points
;; freelimb -- frees the limb darkening parameters in the fits
;; clarlimb -- uses the limb darkening parameters from A. Claret, averaged into wavelength bins
;; psplot -- makes a postscript plot instead of X windows plot
;;           noreject -- No sigma rejection when making plots, shows
;;                       the all the points in the time series
;; differential -- makes a differential measurment the spectrum by
;;                 dividing by one of the bins
;; individual -- plots the individual stars instead of just one at a time
;; pngcopy -- saves an exported PNG file for each PDF plot
;; freeall -- frees the limb darkening parameters, a/R*, impact
;;            parameter as well
;; fixall -- fix all light curve parameters at the values from the
;;           literature (except the DC offset and slope)
;; timebin -- specifies the number of time bins to put the data in

;sigrejcrit = 6D  ;; sigma rejection criterion
sigrejcrit = 5D  ;; sigma rejection criterion
TsigRejCrit = 3D ;; sigma rejection criterion for time bins

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
  endif

  ;; get the binned data
  restore,'data/specdata.sav'

  bingridmiddle = bingrid + binsizes/2E

  nut = n_elements(utgrid)
  nbin = Nwavbins

  ;; get the transit times
  readcol,'transit_info/jan_04_t_time.txt',epoch,tepoch,format='(A,A)',$
          skipline=1

  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
          skipline=1
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor

  u1parm = 0.0E              ;HD 189733
  u2parm = 0.0E

  tstart = date_conv(tepoch[0],'JULIAN')
  tend = date_conv(tepoch[1],'JULIAN')
  tmid = (tend + tstart)/2D

  ;; radius of a planet as a function of wavelength
  plrad = fltarr(nbin)*!values.f_nan
  plrade = fltarr(nbin)*!values.f_nan

  ;; Prepare to save all the planet transit data as a function of wavelength
  paramnames = ['Rp/R*','b_impact','u1','u2','a/R*','linearA','linearB']
  nparams = n_elements(paramnames)
  resultarr = fltarr(nparams,nbin)*!values.f_nan
  resultarrE = resultarr


  ;; orbital phase
  tplot = (utgrid - tmid)/planetdat.period
     
  ;; calculate start and end
  hstart = (tstart - tmid)/planetdat.period
  hend = (tend - tmid)/planetdat.period


  ;; For any binfl error that are zero, set to 0.01
  zerobinp = where(binfle LE 1E-3)
  binfle[zerobinp] = 0.01E

  if n_elements(timebin) NE 0 then begin
     ntime = n_elements(tplot)
     ;; set up the time bins
     timeGrid = (tplot[ntime-1] - tplot[0]) * findgen(timebin)/float(timebin-1l) +$
                tplot[0]
     tsizes = fltarr(timebin) + (tplot[ntime-1l] - tplot[0l])/float(timebin-1l)
     tmiddle = timeGrid + tsizes / 2E
  endif

  for k=0l,nbin-1l do begin
     ;; Reset the x axis (orbital phase, in case it was modified below)
     tplot = (utgrid - tmid)/planetdat.period     
     offp = where(tplot LT hstart OR tplot GT hend)

     if keyword_set(individual) then begin
        y = double(transpose(binind[k,0,*]))
        yerr = double(transpose(binindE[k,0,*]))
        y2 = double(transpose(binind[k,1,*]))
        y2err = double(transpose(binindE[k,1,*]))
        yptitle='Flux (DN)'
        if n_elements(timebin) NE 0 then begin
           print,"Binning not set up for individual star fluxes"
           return
        endif

     endif else begin
        y = double(transpose(binfl[k,*]))
        yerr = double(transpose(binfle[k,*]))
        yptitle='Flux Ratio'
     endelse

     if keyword_set(differential) then begin
;        DiffInd = 9
;        stop
;        y = y / double(transpose(binfl[DiffInd,*]))
        goodwavpts = where((bingridmiddle GT 0.90 and bingridmiddle LT 2.35),ngoodwavpts)
        
        masterY = double(total(binfl[goodwavpts,*],1,/nan))/double(ngoodwavpts)
        y = y/masterY
;        yerr = yerr / double(transpose(binfl[DiffInd,*]))
        yerr = y/masterY
     endif

     meanoff = mean(y[offp],/nan)

     ;; For any binned flux that are NAN, remove
     goodp = where(finite(y) EQ 1)
     if goodp NE [-1] then begin
        y = y[goodp]
        tplot = tplot[goodp]
        offp = where(tplot LT hstart OR tplot GT hend)
     endif
     stdoff = stddev(y[offp])

     rstdoff = robust_sigma(y[offp])
     medoff = median(y[offp])
;     divbycurve = y / 
;     npoly=1
;     result = poly_fit(tplot,y,npoly-1,measure_errors=yerr,yfit=yfit)
;     rgoodp = where(abs(y - medoff) LT rstdoff * Tsigrejcrit,complement=badp)
;     if k EQ 1 then stop


     ;; Throw away points more than n-sigma from the main bunch
     if not keyword_set(noreject) then begin
        for l=0,3-1 do begin    ; iterate 3 times
           stdoff = stddev(y[offp])
           meanoff = mean(y[offp],/nan)
           goodp = where(abs(y - meanoff) LE sigrejcrit * stdoff,complement=throwaframes)
           if goodp NE [-1] then begin
              y = y[goodp]
              tplot = tplot[goodp]
              offp = where(tplot LT hstart OR tplot GT hend)
           endif
        endfor
     endif

     wavname = string(bingridmiddle[k],format='(F4.2)')

     if n_elements(timebin) NE 0 then begin
        ;; Bin the series in time
        ybin = avg_series(tplot,y,y/yerr,timeGrid,tsizes,weighted=1,$
                  oreject=TsigRejCrit,eArr=yerrOut,/silent,errIn=yerr,stdevArr=stdevArr)
        if n_elements(diffInd) NE 0 then begin
           if k EQ diffInd and keyword_set(differential) then begin
              ;; make sure the sigma rejection is off for the
              ;; normalization bin
              ybin = avg_series(tplot,y,y/yerr,timeGrid,tsizes,weighted=1,$
                                eArr=yerrOut,/silent,errIn=yerr,stdevArr=stdevArr)
           endif
        endif

        tplot = tmiddle
        y = ybin
;        yerr = yerrOut
        yerr = stdevArr
        offp = where(tplot LT hstart OR tplot GT hend)
     endif

     if total(finite(y)) GT 0 and total(finite(yerr)) GT 0.0 then begin
        ;; if keyword set, replace the error w/ the off transit stddev
        if keyword_set(offtranserr) then begin
           stdevOff = stddev(y[offp])
           yerr = fltarr(nut) + stdevOff
        endif

        ;find the range where 95% or more of the plots are shown
        if keyword_set(individual) then begin
           ycomb = [y,y2] ;; combine both stars into one array
           sorty = sort(ycomb)
           ylength = n_elements(ycomb)
           ylowerL = ycomb[sorty[ceil(5E/100E*float(ylength))]] * 0.7
           yUpperL = ycomb[sorty[ceil(95E/100E*float(ylength))]] * 1.3
           ydynam = [ylowerL,yUpperL]
        endif else begin
           sorty = sort(y)
           ylength = n_elements(y)
           case 1 of
              keyword_set(fullrange): ydynam = [0,0]
              keyword_set(oneprange): ydynam = [0,1]
              else: begin
                 ylowerL = y[sorty[ceil(5E/100E*float(ylength))]] * 0.97
                 yUpperL = y[sorty[ceil(95E/100E*float(ylength))]] * 1.03
                 ydynam = [ylowerL,yUpperL]
              end
           endcase
        endelse

        if keyword_set(psplot) then begin
           plotnmpre = 'plots/spec_t_series/tser_'+wavname
           device,encapsulated=1, /helvetica,$
                  filename=plotnmpre+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
        endif
;        plot,tplot,y,psym=2,$
        plot,tplot,y,psym=4,$
             xtitle='Orbital Phase',$
             title=wavname+' um Flux ',$
             ytitle=yptitle,$
             yrange=ydynam,ystyle=1
        if keyword_set(individual) then begin
           oplot,tplot,y2,psym=4,color=mycol('blue')
           legend,['Planet Host','Reference Star'],$
                  psym=[2,4],color=mycol(['black','blue']),$
                  /right,/clear
        endif

        ;; print the stdev for y for off points
;        print,'Fractional off transit Stdev in F for ',wavname,': ',stddev(y[offp])/mean(y[offp])
;        print,'Fractional off transit Robust sigma for ',wavname,': ',robust_sigma(y[offp])/median(y[offp])
        ;; try fitting the off transit to a line first
        fitY = linfit(tplot[offp],y[offp])
        Offresid = y[offp]/(fitY[0] + fitY[1]*tplot[offp])
        print,'Frac lin corr robust sigma for ',wavname,': ',robust_sigma(Offresid)/median(Offresid)
        ;; Show the off transit fit
;        oplot,tplot,fity[0] + fity[1]*tplot,color=mycol('red')

        ;; show the transit epochs
        drawy = [!y.crange[0],!y.crange[1]]
        plots,[hstart,hstart],drawy,color=mycol('brown'),linestyle=2
        plots,[hend,hend],drawy,color=mycol('brown'),linestyle=2

        ;; show the off-transit fit for differential measurements
        if keyword_set(differential) then begin
           oplot,tplot,(fitY[0] + fitY[1] *tplot),thick=2,color=mycol('red')
           oploterr,tplot,y,yerr
        endif

        if keyword_set(fitrad) then begin
           ;; fit the data curve
           expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* (P[5] + X * P[6])'
           pi = replicate({fixed:1, limited:[1,0], limits:[0.0E,0.0E]},7)
           pi[0].fixed = 0 ;; make sure the Rp/R* parameter is free
           ;; fix the impact parameter, limb darkening and AoR*
           if keyword_set(freelimb) then begin
              ;; if asked to, free the limb darkening parameter
              pi[2].fixed = 0
              pi[3].fixed = 0
           endif
           if keyword_set(freeall) then begin
              pi[*].fixed = 0
           endif ;; free all parameters
           if keyword_set(fixall) then begin
              pi[*].fixed = 1
           endif
           pi[0].limited = [1,1]
           pi[0].limits = [0D,1D] ;; Keep Rp/R* between 0 and 1
           pi[6].limited = [0,0] ;; let the linear coefficient by + or -
           pi[5].fixed = 0 ;; make sure the flux ratio offset is free
           pi[6].fixed = 0 ;; let the linear coefficient vary
;           if keyword_set(clarlimb) then begin
;              start=[planetdat.p,planetdat.b_impact,u1bin[k],u2bin[k],planetdat.a_o_rstar,0.0E]
;           endif else begin
           start=double([planetdat.p,planetdat.b_impact,u1parm,u2parm,planetdat.a_o_rstar,1.0D,0D])
;           endelse

           result = mpfitexpr(expr,tplot,y,yerr,start,parinfo=pi,perr=punct)
           
           ;; show the fit
           showphase = 0.16E
           nshowpts = 1024      ; number of points to show
           phtest = dindgen(nshowpts)/float(nshowpts)*showphase - showphase/2E
           ytest = quadlc(phtest,result[0],result[1],result[2],result[3],result[4]) $
                   * (result[5] + phtest * result[6])
                                                                                       
           oplot,phtest,ytest,color=mycol('orange'),thick=2

           ;; save the planet radius and all data
           plrad[k] = result[0]
           plrade[k] = punct[0]
           
           resultarr[*,k] = result
           resultarrE[*,k] = punct

           modelY = quadlc(tplot,result[0],result[1],result[2],result[3],result[4]) $
                   * (result[5] + phtest * result[6])

           if keyword_set(psplot) then begin
              device,/close
              device,decomposed=0
              cgPS2PDF,plotnmpre+'.eps',$
                       /delete_ps
              if keyword_set(pngcopy) then begin
                 spawn,'convert -density 160% '+plotnmpre+'.pdf '+plotnmpre+'.png'                 
              endif
              plotnmpre = 'plots/residual_series/residuals_'+wavname
              device,encapsulated=1, /helvetica,$
                     filename=plotnmpre+'.eps'
              device,xsize=14, ysize=10,decomposed=1,/color

           endif

           resid = (y - modelY)/meanoff *100E

           ylowerL = resid[sorty[ceil(5E/100E*float(ylength))]]
           yUpperL = resid[sorty[ceil(95E/100E*float(ylength))]]
           ydynam = [-1E,1E] * max(abs([ylowerL,yUpperL])) * 2E

           plot,tplot,resid,yrange=ydynam,$
                title='Residuals at '+wavname,ystyle=1,$
                xtitle='Orbital Phase',ytitle='Flux Residual (%)',$
                psym=2
;           oploterr,tplot,resid,yerr/meanoff *100E
           drawy = [!y.crange[0],!y.crange[1]]
           plots,[hstart,hstart],drawy,color=mycol('blue'),linestyle=2,thick=2.5
           plots,[hend,hend],drawy,color=mycol('blue'),linestyle=2,thick=2.5
           print,'RMS Residuals (%) for '+wavname,'um',(stddev(y - modelY))/median(y)*100E
;           if wavname EQ '1.14' then stop

        endif
        if keyword_set(psplot) then begin
           device, /close
           device,decomposed=0
           cgPS2PDF,plotnmpre+'.eps',$
                    /delete_ps
           if keyword_set(pngcopy) then begin
              spawn,'convert -density 160% '+plotnmpre+'.pdf '+plotnmpre+'.png'
           endif
        endif
           
;        stop
     endif

  endfor
  
  ;; save the radius data
  
  forprint,bingridmiddle[*],plrad,plrade,$
           textout='radius_vs_wavelength/radius_vs_wavl.txt',$
           comment='#Wavelength  Rp/R*   Rp/R* Error',/silent

  ;; Save the other light curve data
  openw,1,'radius_vs_wavelength/fit_data.txt'
  printf,1,'#Wavelength (um) ',format='(A12,$)'
  for j=0l,nparams-1l do begin
     printf,1,paramnames[j],' Error',format='(A12,A12,$)'
  endfor
  printf,1,''
  for k=0l,nbin-1l do begin
     printf,1,bingridmiddle[k],format='(F12.3,$)'
     for j=0l,nparams-1l do begin
        printf,1,resultarr[j,k],resultarrE[j,k],format='(F12.3,F12.3,$)'
     endfor
     printf,1,''
  endfor
  close,1

  if keyword_set(psplot) then begin
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
end
