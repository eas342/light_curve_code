pro plot_tim_ser,fitcurve=fitcurve,fitpoly=fitpoly,usepoly=usepoly,makestops=makestops,$
                 fullrange=fullrange,smartbin=smartbin,oneprange=oneprange,$
                 offtranserr=offtranserr,freelimbquad=freelimbquad,clarlimb=clarlimb,$
                 psplot=psplot,noreject=noreject,differential=differential,$
                 individual=individual,pngcopy=pngcopy,freeall=freall,fixall=fixall,$
                 timebin=timebin,offreject=offreject,showclipping=showclipping,$
                 errorDistb=errorDistb,colorclip=colorclip,quadfit=quadfit,legorder=legorder,$
                 fixrad=fixrad,freelimblin=freelimblin,showDiffAirmass=showDiffairmass,$
                 nonormalize=nonormalize,showNomRad=showNomRad,fixoffset=fixoffset,$
                 custresidYrange=custresidYrange,fitepoch=fitepoch,singleplot=singleplot,$
                 showmcmc=showmcmc,deletePS=deletePS,showKep=showKep,lindetrend=lindetrend
;; plots the binned data as a time series and can also fit the Rp/R* changes
;; apPlot -- this optional keyword allows one to choose the aperture
;;           to plot
;; fitcurve -- this fitting procedure ueses
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
;; freelimbquad -- frees the quadratic limb darkening parameters in the fits
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
;; offreject -- specifies that the sigma rejection is to be made from
;;              off-transit flux. Otherwise, it uses a model light
;;              curve from literature values
;; showclipping -- shows the clipping of data points to remove outliers
;; errorDistb -- Plot a histogram of the photometric error distribution
;; colorclip -- colors the clipped points
;; quadfit  -- fits a 2nd order Legendre polynomial instead of a linear baseline
;; legorder -- fits a legorder order Legendre polynomial baseline
;; fixrad -- fixes the planet radius to see if limb darkening can do
;;           it all
;; fixoffset -- fixes the out of transit flux at a constant value
;; freelimblin -- frees only the linear limb darkening coefficient
;; nonormalize -- do NOT normalize flux by the off transits points
;; showNomRad -- shows the transit curve with the nominal radius
;;               in addition to the best-fit radius
;; custresidYrange -- sets the Yrange of all residual plots at a
;;                    specific value
;; fitepoch -- fits the transit center
;; singleplot -- Puts everything in a single plot
;; showmcmc -- shows the MCMC results
;; deletePS -- specified whether or not to delete the
;;             postscript file, default is true
;; showKep -- shows the Kepler light curve for KIC 12557548

;sigrejcrit = 6D  ;; sigma rejection criterion
sigrejcrit = 5D  ;; sigma rejection criterion
;TsigRejCrit = 3D ;; sigma rejection criterion for time bins
TsigRejCrit = 2.5D ;; sigma rejection criterion for time bins

if n_elements(deletePS) EQ 0 then deletePS = 1
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
  readcol,'transit_info/transit_epoch.txt',epoch,tepoch,format='(A,A)',$
          skipline=1

  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
          skipline=1
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor

  if keyword_set(showmcmc) then begin
     ;; Read the MCMC parameters
     restore,'data/compiled_model_params.sav' ;; mcmcPars has NxM array 
     ;; N is the wavelength index and M is the parameter index
     ;; get the model evaluation

     modelExpr =''
     openr,1,'data/model_expr.txt'
     readf,1,modelExpr
     close,1
  endif

  u1parm = 0.0E         
  u2parm = 0.0E

  tstart = date_conv(tepoch[0],'JULIAN')
  tend = date_conv(tepoch[1],'JULIAN')
  tmid = date_conv(tepoch[2],'JULIAN')

  ;; radius of a planet as a function of wavelength
  plrad = fltarr(nbin)*!values.f_nan
  plrade = fltarr(nbin)*!values.f_nan

  ;; Prepare to save all the planet transit data as a function of wavelength
;  paramnames = ['Rp/R*','b_impact','u1','u2','a/R*','linearA','linearB','quadC','cubicD']
;  paramnames = ['Rp/R*','b_impact','u1','u2','a/R*','linearA','linearB','quadC','cubicD','quarticE']
  paramnames = ['Rp/R*','b_impact','u1','u2','a/R*','legendre0','legendre1','legendre2','legendre3',$
                'legendre4','legendre5','Phase Offset']
  paramfiles = ['rad'  ,'b_impact','u1','u2','a','legendre0','legendre1','legendre2','legendre3',$
                'legendre4','legendre5','phase_offset']
  nparams = n_elements(paramnames)
  resultarr = fltarr(nparams,nbin)*!values.f_nan
  resultarrE = resultarr

  ;; Prepare to save the off transit RMS
  fracRMSarr = fltarr(nbin)
  fracPhotonarr = fltarr(nbin)

  ;; orbital phase
  tplot = (utgrid - tmid)/planetdat.period
  if keyword_set(showDiffairmass) then begin
     airmassOrig = (1D / sin(!DPI / 180D * altitude[*,0]) - 1D / sin(!DPI / 180D * altitude[*,1]))
  endif else airmassOrig = airmass

  ;; add or subtract integers to phase so that's it's sort of centered
  ;; at 0
  tplot = fold_phase(tplot)

  ;; calculate start and end
  hstart = (tstart - tmid)/planetdat.period
  hend = (tend - tmid)/planetdat.period

  ;; save the phase, time and planet data in array
  save,tmid,tend,tstart,tplot,hstart,hend,$
       planetdat,u1parm,u2parm,$
       filename='data/timedata.sav'

  ;; For any binfl error that are zero, set to 0.01
  zerobinp = where(binfle LE 1E-7)
  if zerobinp NE [-1] then binfle[zerobinp] = 0.01E

  if n_elements(timebin) NE 0 then begin
     ntime = n_elements(tplot)
     ;; set up the time bins
     timeGrid = (tplot[ntime-1] - tplot[0]) * findgen(timebin)/float(timebin) +$
                tplot[0]
     tsizes = fltarr(timebin) + (tplot[ntime-1l] - tplot[0l])/float(timebin)
     tmiddle = timeGrid + tsizes / 2E
  endif

  for k=0l,nbin-1l do begin
     ;; Reset the x axis (orbital phase, in case it was modified below)
     tplot = (utgrid - tmid)/planetdat.period     
     ;; fold back to 0
     tplot = fold_phase(tplot)
     airmass = airmassOrig

     offp = where(tplot LT hstart OR tplot GT hend)

     reffactor=0.45E ;; factor to multiply the reference star by

     if keyword_set(individual) then begin
        y = double(transpose(binind[k,0,*]))
        yerr = double(transpose(binindE[k,0,*]))
        y2 = double(transpose(binind[k,1,*])) * reffactor
        y2err = double(transpose(binindE[k,1,*])) * reffactor
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
        yerr = yerr[goodp]
        if keyword_set(individual) then begin
           y2 = y2[goodp]
           y2err = y2err[goodp]
        endif
        tplot = tplot[goodp]
        airmass = airmass[goodp]
        offp = where(tplot LT hstart OR tplot GT hend)
     endif
     stdoff = stddev(y[offp])

     if not keyword_set(nonormalize) then begin
        y = y / median(y[offp])
        yerr = y / median(y[offp])
     endif

     if not keyword_set(noreject) then begin
        if keyword_set(offreject) then begin
           for l=0,3-1 do begin ; iterate 3 times
              stdoff = stddev(y[offp])
              meanoff = mean(y[offp],/nan)
              goodp = where(abs(y - meanoff) LE sigrejcrit * stdoff,complement=throwaframes)
              if goodp NE [-1] then begin
                 y = y[goodp]
                 yerr = yerr[goodp]
                 tplot = tplot[goodp]
                 if keyword_set(individual) then begin
                    y2 = y2[goodp]
                    y2err = y2err[goodp]
                 endif
                 airmass = airmass[goodp]
                 offp = where(tplot LT hstart OR tplot GT hend)
              endif
           endfor
        endif else begin
           if keyword_set(differential) then divbycurve = y else begin
              ;; before doing the sigma rejection, divide by literature model
              ;; for 
              ymodel = quadlc(tplot,planetdat.p,planetdat.b_impact,$
                              u1parm,u2parm,planetdat.a_o_rstar)
              divbycurve = y / ymodel
           endelse
           tplotdivcurves = tplot ;; if plotting the pre-rejection data, save the independent variable

           ;throw away all n_sigma events before de-trending
           firstCutSig = 12E
           rstdoff = robust_sigma(y[offp])
           medoff = median(y[offp])

           goodp = where(abs(y - medoff) LE firstCutSig * rstdoff,complement=throwaways)
           if goodp NE [-1] then begin
              if throwaways NE [-1] then begin
                 tclip1 = tplot[throwaways]
                 yclip1 = y[throwaways]
              endif
              yfull = y
              y = y[goodp]
              divbycurveclip1 = divbycurve[goodp]
              yerr = yerr[goodp]
              if keyword_set(individual) then begin
                 y2 = y2[goodp]
                 y2err = y2err[goodp]
              endif
              tplot = tplot[goodp]
              airmass = airmass[goodp]
           endif
           ;; fit result to a robust line
           rlinefit = robust_linefit(tplot,divbycurveclip1,yfit)
           ;; divide by the line to flatten out
           yflat = divbycurveclip1 / yfit
           rsigma = robust_sigma(yflat)
           maxval = max(yflat,maxp)
           minval = min(yflat,minp)
           nombysigma = (yflat - mean(yflat))/rsigma
           if keyword_set(lindetrend) then y = y/yfit
;           secondCutsig = 4.0E
           secondCutsig = 3.5E

           goodp = where(abs(nombysigma) LE secondCutsig,complement=throwaways2)
           if goodp NE [-1] then begin
              if throwaways2 NE [-1] then begin
                 yclip2 = y[throwaways2]
                 tclip2 = tplot[throwaways2]
              endif
              y = y[goodp]
;              yerr = fltarr(n_elements(goodp)) + rsigma * rlinefit[0]
              yerr = yerr[goodp]
              if keyword_set(individual) then begin
                 y2 = y2[goodp]
                 y2err = y2err[goodp]
              endif
              tplot = tplot[goodp]
              airmass = airmass[goodp]
              offp = where(tplot LT hstart OR tplot GT hend)
           endif

        endelse


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

        airbin = avg_series(tplot,airmass,fltarr(n_elements(airmass)),timeGrid,tsizes,$
                            weighted=0,/silent)
        tplot = tmiddle

        y = ybin
        fracPhotonarr[k] = yerrOut[0] / ybin[0]
;        yerr = yerrOut
        yerr = stdevArr
        airmass = airbin
        offp = where(tplot LT hstart OR tplot GT hend)
        if keyword_set(offtranserr) then begin
           fitY2 = linfit(tplot[offp],y[offp])
           Offresid2 = y[offp] - (fitY2[0] + fitY2[1]*tplot[offp])
           rstdevOff2 = robust_sigma(Offresid2)
           yerr = fltarr(n_elements(yerr)) + rstdevOff2
        endif

     endif

     if total(finite(y)) GT 0 and total(finite(yerr)) GT 0.0 then begin
        ;; if keyword set, replace the error w/ the off transit stddev

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
              keyword_set(differential): begin
                 ylowerL = y[sorty[ceil(5E/100E*float(ylength))]] * 0.999
                 yUpperL = y[sorty[floor(95E/100E*float(ylength))]] * 1.002
                 ydynam = [ylowerL,yUpperL]
              end
              else: begin
                 ylowerL = y[sorty[ceil(5E/100E*float(ylength))]] * 0.99
                 yUpperL = y[sorty[floor(95E/100E*float(ylength))]] * 1.01
                 ydynam = [ylowerL,yUpperL]
              end
           endcase
        endelse

        if keyword_set(psplot) and (not keyword_set(singleplot) OR k EQ 0)  then begin
           plotnmpre = 'plots/spec_t_series/tser_'+wavname
           device,encapsulated=1, /helvetica,$
                  filename=plotnmpre+'.eps'
           if keyword_set(singleplot) then begin
              device,xsize=11, ysize=14.1,decomposed=1,/color
           endif else device,xsize=14, ysize=10,decomposed=1,/color
        endif
;        plot,tplot,y,psym=2,$
        custXrange=[-0.1,0.1]
        if keyword_set(singleplot) then begin
           if k EQ 0 then begin
              ;; Set up everything for the first time plot
              myNoerase=0
              yptitle= yptitle + ' + Offset'
              myXtitle='Orbital Phase'
              tickformat='(G0)'
              myXrange=[min(tplot),max(tplot)+0.25*(max(tplot)-min(tplot))]
              ;; If it's Dec 23, make an adjustment to avoid
              ;; squishing numbers
              if myXrange[0] GT -0.06254 and myXrange[0] LE -0.0625 and $
                 myXrange[1] GT 0.0954 and myXrange[1] LE 0.09542 then myXrange=[-0.07,0.1]
           endif else begin
              ;; Set up for subsequent plots
              myNoerase=1
              yptitle=''
              myXtitle=''
              tickformat='(A1)'
           endelse
           spacing=0.025
           ydynam=[1E - spacing * (1+Nwavbins),1+spacing]
           offset = k * spacing
           myTitle=''
        endif else begin
           myNoErase=0
           offset = 0
           myTitle=wavname+' um Flux'
           tickformat='(G0)'
           myXtitle='Orbital Phase'
        endelse
        plot,tplot,y,psym=4,$
             xtitle=myXtitle,$
             title=myTitle,$
             ytitle=yptitle,$
             yrange=ydynam,ystyle=1,/nodata,xstyle=1,$
             noerase=myNoErase,$
             xtickformat=tickformat,ytickformat=tickformat,$
             xrange=myXrange
        if not keyword_set(differential) then begin
           if keyword_set(showclipping) then begin
              oplot,tplot,y,psym=5,color=mycol('red') ;; original data
              oplot,tplotdivcurves,divbycurve,psym=4  ;; divided by light curve
              oplot,tplot,yfit,color=mycol('blue')    ;; fitted line to curve
           endif else begin
              if n_elements(timebin) EQ 0 then oplot,tplot,y-offset,psym=4 else begin
                 oploterror,tplot,y-offset,tsizes/2E,yerr,psym=3,hatlength=0,thick=2
              endelse
           endelse 
        endif
        if keyword_set(individual) then begin
           oplot,tplot,y2,psym=4,color=mycol('blue')
           legend,['Planet Host','Reference Star X '+strtrim(reffactor,1)],$
                  psym=[4,4],color=mycol(['black','blue']),$
                  /right,/clear
        endif
        if keyword_set(singleplot) then begin
           xyouts,!x.crange[1]-0.1*(!x.crange[1]-!x.crange[0]),$
                  median(y)-offset,$
                  [wavname+' um'],alignment=0.5
        endif
        ;; print the stdev for y for off points
;        print,'Fractional off transit Stdev in F for ',wavname,': ',stddev(y[offp])/mean(y[offp])
;        print,'Fractional off transit Robust sigma for ',wavname,': ',robust_sigma(y[offp])/median(y[offp])
        ;; try fitting the off transit to a function first
        if keyword_set(quadfit) then begin
           result = poly_fit(tplot[offp],y[offp],2,measure_errors=yerr[offp],yfit=yfit)
           Offresid = y[offp] - yfit
        endif else begin
           fitY = linfit(tplot[offp],y[offp])
           Offresid = y[offp] - (fitY[0] + fitY[1]*tplot[offp])
        endelse
        fracRMSarr[k] = robust_sigma(Offresid)/median(y[offp])
        print,'Frac lin corr robust sigma for ',wavname,': ',fracRMSarr[k]
        ;; Show the off transit fit
;        oplot,tplot,fity[0] + fity[1]*tplot,color=mycol('red')
        if keyword_set(offtranserr) then begin
           rstdevOff = robust_sigma(Offresid)
           yerr = fltarr(n_elements(goodp)) + rstdevOff
        endif

        ;; show the transit epochs
        drawy = [!y.crange[0],!y.crange[1]]
        plots,[hstart,hstart],drawy,color=mycol('brown'),linestyle=2
        plots,[hend,hend],drawy,color=mycol('brown'),linestyle=2

        ;; show the off-transit fit for differential measurements (as
        if keyword_set(differential) then begin
           if not keyword_set(fitcurve) then begin
                 ;; long as fits aren't being performed)
                 
              if keyword_set(quadfit) then begin
                 oplot,tplot,(result[0] + result[1] *tplot + result[2] * tplot^2),$
                       thick=2,color=mycol('red')
              endif else begin
                 oplot,tplot,(fitY[0] + fitY[1] *tplot),thick=2,color=mycol('red')
              endelse
           endif
           if n_elements(timebin) EQ 0 then tsizes = fltarr(n_elements(tplot)) + tplot[1]-tplot[0]
           oploterror,tplot,y,tsizes/2E,yerr,psym=3,hatlength=0,thick=2
        endif
        ;;plot the clipped points
        if keyword_set(colorclip) then begin
           if throwaways NE [-1] then oplot,tclip1,yclip1,psym=6,color=mycol('blue')
           if throwaways2 NE [-1] then oplot,tclip2,yclip2,psym=5,color=mycol('blue')
        endif
        if keyword_set(makestops) then stop
        
        if keyword_set(errorDistb) then begin
           if keyword_set(psplot) then begin
              device,/close
              device,decomposed=0
              cgPS2PDF,plotnmpre+'.eps',$
                       delete_ps=deletePS
              if keyword_set(pngcopy) then begin
                 spawn,'convert -density 160% '+plotnmpre+'.pdf '+plotnmpre+'.png'
              endif
              plotnmpre = 'plots/error_distrib/error_hist_'+wavname
              device,encapsulated=1, /helvetica,$
                     filename=plotnmpre+'.eps'
              device,xsize=14, ysize=10,decomposed=1,/color

           endif
           binsize=0.25E
           plothist,nombysigma,xhist,yhist,bin=binsize,$ ;; show error distribution
                    ytitle='Number of Points',$
                    xtitle='(Flux - median)/Robust Sigma'
           xgaussian = findgen(100)/(100E - 1E) * (!x.crange[1] - !x.crange[0])+$
                       !x.crange[0]
           totY = total(yhist)*binsize
           ygaussian = gaussian(xgaussian,[totY/sqrt(2E * !DPI),0,1])
           oplot,xgaussian,ygaussian,color=mycol('red')

        endif

        if keyword_set(showmcmc) then begin
           mcmcShowP = 350
           phaseShow = findgen(mcmcshowP)/float(mcmcShowP) * (max(tplot)-min(tplot)) + min(tplot)
           ;; First find the mean function with expression evaluation
           meanfunctest = expression_eval(modelExpr,phaseShow,mcmcPars[k,0:8])
           ;oplot,phaseShow,meanfunctest-offset,color=mycol('blue'),thick=2
           ;; Find the residual vector
           meanfuncdat = expression_eval(modelExpr,float(tplot),mcmcPars[k,0:8])
           mcmcResid = y - meanfuncdat
           ;; Find the inverse covariance matrix using the likelihood
           ;; function which also needs the inverse covariance matrix
           likeL = ev_leval([mcmcPars[k,9],mcmcPars[k,10],0],x=tplot,yin=y,yerr=yerr,Cinv=Cinv)
           mcmcModel = fltarr(n_elements(phaseShow))
           for m=0l,mcmcShowP-1l do begin
              columnvec = fltarr(n_elements(tplot))
              Argument = -abs((phaseShow[m] - tplot))*mcmcPars[k,10]
              evalpmcmc = where(Argument GT -15D) ;; ensure the argument of exponential doesn't cause overflow error
              if evalpmcmc NE [-1] then $
;                 columnvec[evalpmcmc] = mcmcPars[k,9]^2 * exp(Argument[evalpmcmc]) * (1E + abs(phaseShow[m] -  tplot) * mcmcPars[k,10])
                 columnvec[evalpmcmc] = mcmcPars[k,9] * exp(Argument[evalpmcmc]) / mcmcPars[k,10]
              mcmcModel[m] = meanfunctest[m] + columnvec ## Cinv ## transpose(mcmcResid)
           endfor
           oplot,phaseShow,mcmcModel-offset,color=mycol('blue'),thick=2
        endif

        if keyword_set(showKep) then begin
           ;; Show the Kepler Light Curve
           readcol,'data/phase_folded_kepler_deep_kic1255.csv',kphaseD,kfluxD,skipline=6,$
                   format='(F,F)',/silent
           readcol,'data/phase_folded_kepler_shallow_kic1255.csv',kphaseS,kfluxS,skipline=6,$
                   format='(F,F)',/silent
           oplot,kphaseD,kfluxD-offset,color=mycol('red')
           oplot,kphaseS,kfluxS-offset,color=mycol('red')

        endif

        if keyword_set(fitcurve) then begin
           ;; fit the data curve

;           start=double([planetdat.p,planetdat.b_impact,u1parm,u2parm,$
;                         planetdat.a_o_rstar,1.0D,0D,0D,0D])
           start=double([planetdat.p,planetdat.b_impact,u1parm,u2parm,$
                         planetdat.a_o_rstar,1.0D,0D,0D,0D,0D,0D,0D])


;           if keyword_set(quadfit) then begin
           if keyword_set(differential) then begin
              quadlcArg ='X'
              for m=0l,4 do begin
                 quadlcArg=quadlcArg+','+strtrim(start[m],1)
              endfor
              expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* (P[5] + X * P[6] + X^2 * P[7] + X^3'+$
                     ' * P[8])/quadlc('+quadlcArg+')'

           endif else begin
              expr = 'quadlc(X-P[11],P[0],P[1],P[2],P[3],P[4])* ( P[5] + '+$
                     'Legendre((2D * X - Max(X) - Min(X))/(Max(X) - Min(X) + 3D-16),1) * P[6] + '+$
                     'Legendre((2D * X - Max(X) - Min(X))/(Max(X) - Min(X) + 3D-16),2) * P[7] + '+$
                     'Legendre((2D * X - Max(X) - Min(X))/(Max(X) - Min(X) + 3D-16),3) * P[8] + '+$
                     'Legendre((2D * X - Max(X) - Min(X))/(Max(X) - Min(X) + 3D-16),4) * P[9] + '+$
                     'Legendre((2D * X - Max(X) - Min(X))/(Max(X) - Min(X) + 3D-16),5) * P[10])'
              
              ;expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* (P[5] + X * P[6] + X^2 * P[7] + X^3 * P[8])'
           endelse
;              pi = replicate({fixed:1, limited:[1,0], limits:[0.0E,0.0E]},8)
;           endif else begin
;              expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* (P[5] + X * P[6])'
;           pi = replicate({fixed:1, limited:[1,0],limits:[0.0E,0.0E]},9)
           pi = replicate({fixed:1, limited:[1,0], limits:[0.0E,0.0E]},nparams)
;           endelse
           ;; make sure the Rp/R* parameter is free
           if not keyword_set(fixrad) then pi[0].fixed = 0 
           ;; fix the impact parameter, limb darkening and AoR*
           case 1 of
              keyword_set(freelimbquad): begin
                 ;; if asked to, free the quadratic limb darkening parameter
                 pi[2].fixed = 0
                 pi[3].fixed = 0
              end
              keyword_set(freelimblin): pi[2].fixed = 0
              else: junk=junk
           endcase

           if keyword_set(freeall) then begin
              pi[*].fixed = 0
           endif ;; free all parameters
           if keyword_set(fixall) then begin
              pi[*].fixed = 1
           endif
           pi[0].limited = [1,1] ;; make sure Rp/R* is limited
           pi[0].limits = [0D,1D] ;; Keep Rp/R* between 0 and 1
           ;; make sure the flux ratio offset is free
           if keyword_set(fitepoch) then begin
              pi[11].fixed = 0
              pi[11].limited = [0,0]
           endif
           if not keyword_set(fixoffset) then pi[5].fixed = 0 
           ;; Let the limb darkening be + or -
           pi[2].limited = [0,0]
           pi[3].limited = [0,0]
           for i=6,10-1 do pi[i].limited = [0,0] ;; let the polynomial coefficients be + or -

           case 1 of 
              keyword_set(quadfit): begin
                 pi[6].fixed = 0 ;; let the linear coefficient
                 pi[7].fixed = 0 ;; let the second Legendre coefficient vary
              end
              n_elements(legOrder) NE 0: begin
                 for i=5l,legOrder+5l do begin
                    pi[i].fixed = 0 ;; let the i-th Legendre coefficient vary
                 endfor
              end
              else: pi[6].fixed = 0 ;; let the linear coefficient vary by default
           endcase

;           if keyword_set(clarlimb) then begin
;              start=[planetdat.p,planetdat.b_impact,u1bin[k],u2bin[k],planetdat.a_o_rstar,0.0E]
;           endif else begin

;           endelse
;           stop
           result = mpfitexpr(expr,tplot,y,yerr,start,parinfo=pi,perr=punct)
           modelPts = 512l
           ntpoints = n_elements(tplot)
           modelY = expression_eval(expr,tplot,result)
           modelX = findgen(modelPts)/(modelPts-1l) * (tplot[ntpoints-1] - tplot[0]) + tplot[0]
           modelY1 = expression_eval(expr,modelX,result)
           oplot,modelX,modelY1-offset,color=mycol('blue'),thick=2

           ;; save the planet radius and all data
           plrad[k] = result[0]
           plrade[k] = punct[0]
           
           resultarr[*,k] = result
           resultarrE[*,k] = punct
           
           if keyword_set(showNomRad) then begin
              pi[0].fixed = 1 ;; fix the radius
              nomRadResult = mpfitexpr(expr,tplot,y,yerr,start,parinfo=pi)
              modelY2 = expression_eval(expr,modelX,nomRadResult)
              oplot,modelX,modelY2,color=mycol('orange'),thick=2
              legend,['Best-Fit Radius','Nominal Radius'],$
                     linestyle=[0,0],color=mycol(['blue','orange']),$
                     thick=[2,2],/clear
           endif

           resid = (y - modelY)/meanoff *100E

           if not keyword_set(singleplot) then begin
              if keyword_set(psplot) then begin
                 device,/close
                 device,decomposed=0
                 cgPS2PDF,plotnmpre+'.eps',$
                          delete_ps=deletePS
                 if keyword_set(pngcopy) then begin
                    spawn,'convert -density 160% '+plotnmpre+'.pdf '+plotnmpre+'.png'                 
                 endif
                 plotnmpre = 'plots/residual_series/residuals_'+wavname
                 device,encapsulated=1, /helvetica,$
                        filename=plotnmpre+'.eps'
                 device,xsize=14, ysize=10,decomposed=1,/color
                 
              endif
              
              ylowerL = resid[sorty[ceil(5E/100E*float(ylength))]]
              yUpperL = resid[sorty[floor(95E/100E*float(ylength))]]
              ydynam = [-1E,1E] * max(abs([ylowerL,yUpperL])) * 4E
              
              if n_elements(custresidYrange) NE 0 then ydynam = custresidYrange
              
              overplotMarg = [13,14]
              plot,tplot,resid,yrange=ydynam,$
                   title='Residuals at '+wavname,$
                   xtitle='Orbital Phase',ytitle='Flux Residual (%)',$
                   psym=2,ystyle=8+1,xmargin=overplotMarg,/nodata,$
                   xrange=custXrange
              if n_elements(timebin) EQ 0 then oplot,tplot,resid,psym=4 else begin
                 oploterror,tplot,resid,tsizes/2E,yerr/meanoff * 100E,psym=3,hatlength=0,thick=2
              endelse
              
              prevXrange=!x.crange
              plot,tplot,airmass,xstyle=1+4,xrange=prevXrange,$
                   /noerase,ystyle=4+16,/nodata,xmargin=overplotMarg
              oplot,tplot,airmass,color=mycol('blue')
              if keyword_set(showDiffairmass) then begin
                 airmassname = 'Differential Airmass'
              endif else airmassname = 'Airmass'
              axis,yaxis=1,yrange=!y.crange,ystyle=1,color=mycol('blue'),$
                   ytitle=airmassname
;           oploterr,tplot,resid,yerr/meanoff *100E
              drawy = [!y.crange[0],!y.crange[1]]
              plots,[hstart,hstart],drawy,color=mycol('blue'),linestyle=2,thick=2.5
              plots,[hend,hend],drawy,color=mycol('blue'),linestyle=2,thick=2.5
              print,'RMS Residuals (%) for '+wavname,'um',(stddev(y - modelY))/median(y)*100E
           endif
        endif
        if not keyword_set(fitcurve) then begin
           ntplot = n_elements(tplot)
           modelY = fltarr(ntplot)
           resid = fltarr(ntplot)
        endif
        forprint,tplot,y,yerr,modelY,resid,$
                 textout='data/cleaned_tim_ser/timeser_'+wavname+'um_.txt',$
                 comment='#Phase  Flux  Fl_err  Model_fl   Residual for '+wavname+'um',$
                 /silent

        if keyword_set(psplot) and (not keyword_set(singleplot) OR k EQ Nwavbins-1l)  then begin
           device, /close
           device,decomposed=0
           cgPS2PDF,plotnmpre+'.eps',$
                    delete_ps=deletePS
           if keyword_set(pngcopy) then begin
              if keyword_set(singleplot) then begin
                 spawn,'convert -density 300% '+plotnmpre+'.pdf '+plotnmpre+'.png'
              endif else begin
                 spawn,'convert -density 160% '+plotnmpre+'.pdf '+plotnmpre+'.png'
              endelse
           endif
        endif

     endif

  endfor
  
  ;; save the radius data
  
  forprint,bingridmiddle[*],binsizes,plrad,plrade,$
           textout='radius_vs_wavelength/radius_vs_wavl.txt',$
           comment='#Wavelength(um) Binsize (um)  Rp/R*   Rp/R* Error',/silent

  ;; save the RMS of off transit fluxu
  save,bingrid,fracRMSarr,tsizes,bingridmiddle,binsizes,$
       planetdat,fracPhotonarr,$
       filename='data/rmsdata.sav'

  ;; Save the parameters to file
  for j=0l,nparams-1l do begin
     forprint,bingridmiddle[*],binsizes,resultarr[j,*],resultarrE[j,*],$
              textout='radius_vs_wavelength/fit_data/'+paramfiles[j]+'_vs_wavl.txt',$
              comment=string('#Wavelength(um)','Binsize (um)',paramnames[j],'Error',$
                             format='(4(A16))'),$
              format='(4(G16.8))',/silent
     
  endfor
  openw,1,'radius_vs_wavelength/fit_data.txt'
  printf,1,'#Wavelength (um) ',format='(A16,$)'
  printf,1,'Bin Size (um)',format='(A16,$)'
  for j=0l,nparams-1l do begin
     printf,1,paramnames[j],' Error',format='(A16,A16,$)'
  endfor
  printf,1,''
  for k=0l,nbin-1l do begin
     printf,1,bingridmiddle[k],format='(F16.8,$)'
     printf,1,binsizes[k],format='(F16.8,$)'
     for j=0l,nparams-1l do begin
        printf,1,resultarr[j,k],resultarrE[j,k],format='(F16.8,F16.8,$)'
     endfor
     printf,1,''
  endfor
  close,1

  if keyword_set(singleplot) then begin
     !p.multi = 0
     !Y.Omargin = [0,0]
  endif

  if keyword_set(psplot) then begin
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
end
