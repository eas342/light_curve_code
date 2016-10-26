pro plot_tim_ser,fitcurve=fitcurve,fitpoly=fitpoly,usepoly=usepoly,makestops=makestops,$
                 fullrange=fullrange,smartbin=smartbin,oneprange=oneprange,$
                 offtranserr=offtranserr,straightupRMS=straightupRMS,$
                 freelimbquad=freelimbquad,clarlimb=clarlimb,$
                 psplot=psplot,noreject=noreject,differential=differential,$
                 individual=individual,background=background,ratioback=ratioback,$
                 pngcopy=pngcopy,freeall=freeall,fixall=fixall,$
                 timebin=timebin,offreject=offreject,showclipping=showclipping,$
                 errorDistb=errorDistb,colorclip=colorclip,quadfit=quadfit,legorder=legorder,$
                 fixrad=fixrad,freelimblin=freelimblin,showDiffAirmass=showDiffairmass,$
                 nonormalize=nonormalize,showNomRad=showNomRad,fixoffset=fixoffset,$
                 custresidYrange=custresidYrange,fitepoch=fitepoch,singleplot=singleplot,$
                 showmcmc=showmcmc,deletePS=deletePS,showKep=showKep,lindetrend=lindetrend,$
                 showjump=showjump,kepfit=kepfit,skipReset=skipReset,custSep=custSep,$
                 showNomMCMC=showNomMCMC,useGPasfit=useGPasfit,kepdiff=kepdiff,$
                 custyrange=custyrange,tryAlt=tryAlt,trycorrect=trycorrect,$
                 secondary=secondary,$
                 presentation=presentation,slitmod=slitmod,fixprof=fixprof,psmooth=psmooth,$
                 custxrange=custxrange,noplots=noplots,custTitle=custTitle,tmedian=tmedian,$
                 custxmargin=custxmargin,custymargin=custymargin,labelKep=labelKep,boot=boot,$
                 skipwavl=skipwavl,sinfit=sinfit,showBJD=showBJD,wavelabelcsize=wavelabelcsize
;; plots the binned data as a time series and can also fit the Rp/R* changes
;; apPlot -- this optional keyword allows one to choose the aperture
;;           to plot
;; fitcurve -- this fitting procedure ueses
;;            planet transit information to fit Rp/R* as a
;;            function of wavelength
;; fitkep -- fit the KIC 1255 light curve using Kepler curves
;; fitpoly -- this fits a polynomial to the data to take out long term
;;           trends
;; usepoly -- this uses the polynomial fit from other sources
;;            (calibrators) on a given source
;; npoly -- allows one to set the number of polynomial parameters used
;; fullrange -- allows for plotting of the full y range
;; smartbin -- use the smart-binned data instead of the binned data
;; oneprange -- sets all the y ranges to be 1%
;; offtranserr -- makes the error equal to the stddev of
;;                linearly-detrended off-transit points
;; straightupRMS -- with offtranserr set, it uses the standard
;;                  deviation of out of transit flux (but no linear
;;                  de-trending applied)
;; freelimbquad -- frees the quadratic limb darkening parameters in the fits
;; clarlimb -- uses the limb darkening parameters from A. Claret, averaged into wavelength bins
;; psplot -- makes a postscript plot instead of X windows plot
;;           noreject -- No sigma rejection when making plots, shows
;;                       the all the points in the time series
;; differential -- makes a differential measurment the spectrum by
;;                 dividing by one of the bins
;; individual -- plots the individual stars instead of just one at a
;;               time
;; individual -- plots the stars' backgrounds instead of just one at a time
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
;; showNomMCMC - show the MCMC fit without GP hyperparameters
;; deletePS -- specified whether or not to delete the
;;             postscript file, default is true
;; showKep -- shows the Kepler light curve for KIC 12557548
;; showjump -- show where the telescope tracking jump occured
;; skipReset - skip the reset at the end of the program
;; custSep -- custom separation between data in single plot mode
;; useGPasfit -- use the Gaussian process model as a fit (instead of
;;               just the mean function) The result should be like
;;               white noise
;; kepfit -- use the KIC 12557458b light curve fit (scaled from the
;;           Kepler Phot)
;; kepdiff -- same as kepfit but for differential light curve fitting
;; custyrange -- custom y range for time series plots
;; tryAlt - try an alternative time series that has been corrected to
;;          a non-linear trend between stars 1 and 2
;; trycorrect - try a corrected time series where a model using the state
;;              parameters was fit to the data and then divided out
;; secondary -- designed for secondary eclipse
;; presentation -- makes things bigger for a power point presentation
;; slitmod - Use a slit loss model
;; fixpro - fix the profile of the stars in the slit loss model
;; psmooth - sets the smooth size used in smothing the optical state
;;           parameters in the slit model
;; noplots - skips all plots (just collects data/fits)
;; custTitle - allows you to specify a custom title
;; tmedian - do median combining instead of average weighting
;; labelKep - label the Kepler average light curve
;; boot - find bootstrap errors instead of mpfit
;; skipwavl - skip the wavelength label
;; showBJD - show the BJD at the top X axis
;; wavelabelcsize - custom wavelength label character size

;sigrejcrit = 6D  ;; sigma rejection criterion
sigrejcrit = 5D  ;; sigma rejection criterion
;TsigRejCrit = 3D ;; sigma rejection criterion for time bins
TsigRejCrit = 5D ;; sigma rejection criterion for time bins

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
  readcol,'transit_info/transit_epoch.txt',epoch,tepoch,format='(A,D)',$
          skipline=1,/silent

  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
          skipline=1,/silent
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor

  ;; Get the (plot spacing) separation between different time series for a given planet
  readcol,'param_input/time_series_sep.txt',separationA,format='(F)',skipline=1,/silent


  if keyword_set(showmcmc) then begin
     ;; Read the MCMC parameters
     restore,'data/compiled_model_params.sav' ;; mcmcPars has NxM array 
     ;; N is the wavelength index and M is the parameter index
     ;; get the model evaluation


     modelExpr =''
     openr,1,'data/model_expr.txt'
     readf,1,modelExpr
     close,1
     readcol,'param_input/kernel_choices.txt',$
             skipline=1,format='(A,A)',$
             instrumentkernRef,kernchoice
  endif

  if keyword_set(presentation) then begin
     PSplotXSize = 9
     PSplotYSize = 6
     PSSingleXsize = 9
     PSSingleYsize = 6
  endif else begin
     PSplotXSize = 14
     PSplotYSize = 10
     PSSingleXsize = 12
     PSSingleYsize = 9
  endelse

  u1parm = 0.0E         
  u2parm = 0.0E

  tstart = tepoch[0]
  tend = tepoch[1]
  tmid = tepoch[2]

  ;; radius of a planet as a function of wavelength
  plrad = fltarr(nbin)*!values.f_nan
  plrade = fltarr(nbin)*!values.f_nan

  ;; Prepare to save all the planet transit data as a function of
  ;; wavelength
  case 1 of
     keyword_set(slitmod): begin
        paramnames = ['Rp/R*','b_impact','u1','u2','a/R*','c_0','c_2',$
                      'B','H','Phase Offset','legendre0','legendre1','legendre2']
        paramfiles = ['rad'  ,'b_impact','u1','u2','a',   'c_0','c_2',$
                      'B','H','phase_offset','legendre0','legendre1','legendre2']
     end
     keyword_set(sinfit): begin
        paramnames = ['Amp','period','phaseoffset','blank03','blank04',$
                      'blank05','blank06','blank07','legendre0','legendre1','legendre2']
        paramfiles = ['mmp','period','phase_offset','blank03','blank04',$
                      'blank05','blank06','blank07','legendre0','legendre1','legendre2']
     end
     else: begin
        paramnames = ['Rp/R*','b_impact','u1','u2','a/R*','legendre0','legendre1','legendre2','legendre3',$
                      'legendre4','legendre5','Phase Offset']
        paramfiles = ['rad'  ,'b_impact','u1','u2','a','legendre0','legendre1','legendre2','legendre3',$
                      'legendre4','legendre5','phase_offset']
     end
  endcase
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
  tplot = fold_phase(tplot,secondary=secondary)

  if keyword_set(custxrange) then begin
     lookP = where(tplot GT custxrange[0] and tplot LT custxrange[1])
     tplot = tplot[lookp]
     binfl = binfl[*,lookp]
     binfle = binfle[*,lookp]
     if n_elements(binind) NE 0 then begin
        binind = binind[*,*,lookp]
        binindE = binindE[*,*,lookp]
     endif
     airmass = airmass[lookp]
     utgrid = utgrid[lookp]
  endif

  ;; calculate start and end
  hstart = (tstart - tmid)/planetdat.period
  hend = (tend - tmid)/planetdat.period

  restore,'data/used_date.sav'
  tdataUseDate = usedate ;; save the speclist used date associated with the file
  ;; save the phase, time and planet data in array
  save,tmid,tend,tstart,tplot,hstart,hend,$
       planetdat,u1parm,u2parm,tdataUseDate,$
       filename='data/timedata.sav'

  ;; try the alternate flux grid
  if keyword_set(tryAlt) then begin
     two_star_corr
     restore,'data/alt_tim_ser.sav'
     binfl = fltarr(1,n_elements(tplot))
     binfl[0,*] = altY
  endif

  if keyword_set(trycorrect) then begin
     restore,'data/state_parameters/alt_tim_ser.sav'
     binfl = fltarr(1,n_elements(tplot))
     binfl[0,*] = yfullcorrected
  endif

  if keyword_set(slitmod) then begin
     ;; Get the spec list Prefix name
     restore,'data/used_date.sav'
     ;; Get the state parameters
     restore,'data/state_parameters/full_parameters/'+$
             specfileListNamePrefix+'.sav'
     npts = n_elements(tplot)
     nstate = 7 ;; number of state parameters
     inputX = dblarr(nstate,npts)
     inputX[0,*] = statePstruct.specshift1
     inputX[1,*] = statePstruct.specshift2
     inputX[2,*] = statePstruct.fwhm1
     inputX[3,*] = statePstruct.fwhm2
     inputX[4,*] = statePstruct.phase
     inputX[5,*] = statePstruct.voigtdamp1
     inputX[6,*] = statePstruct.voigtdamp2
     star1ind = [2,5]
     star2ind = [3,6]
     nstarts = n_elements(star1ind)
     for i=0l,nstarts-1l do begin
        avgX = (inputX[star1ind[i],*] + inputX[star2ind[i],*])/2E
        inputX[star1ind[i],*] = avgX
        inputX[star2ind[i],*] = avgX
        if keyword_set(fixProf) then begin
           inputX[star1ind[i],*] = median(inputX[star1ind[i],*])
           inputX[star2ind[i],*] = median(inputX[star2ind[i],*])
        endif
     endfor
     ;; Prepare to reject all n sigma outliers in state parameters
     StatesigthreshS=5E
     badArray = intarr(npts)
     for i=0l,nstate-1 do begin
        if i NE 4l then begin
           tempArray = inputX[i,*]
           rsigma = robust_sigma(tempArray)
           goodp = where(abs(tempArray - median(tempArray)) LT StatesigthreshS * rsigma,complement=badp)
           if badP NE [-1] then begin
              badArray[badp] = 1
           endif
        endif
     endfor
     inputXorig = inputX
  endif else begin
     badArray = intarr(n_elements(tplot)) ;; no bad points
  endelse

  ;; For any binfl error that are zero, set to 0.01
  zerobinp = where(binfle LE 1E-7)
  if zerobinp NE [-1] then binfle[zerobinp] = 0.01E

  if n_elements(timebin) NE 0 then begin
     ntime = n_elements(tplot)
     ;; set up the time bins
     timeGrid = (max(tplot) - min(tplot)) * findgen(timebin)/float(timebin) +$
                min(tplot)
     tsizes = fltarr(timebin) + (max(tplot) - min(tplot))/float(timebin)
     tmiddle = timeGrid + tsizes / 2E
  endif

  for k=0l,nbin-1l do begin
     ;; Reset the x axis (orbital phase, in case it was modified below)
     tplot = (utgrid - tmid)/planetdat.period     
     ;; fold back to 0
     tplot = fold_phase(tplot,secondary=secondary)
     airmass = airmassOrig
     ;; Same reset applies to inputX
     if keyword_set(slitmod) then begin
        inputX = inputXorig
     endif
     offp = where(tplot LT hstart OR tplot GT hend)

     reffactor=0.45E ;; factor to multiply the reference star by

     case 1 of
        keyword_set(ratioback): begin
           y = double(transpose(binback[k,0,*]))
           yerr = double(transpose(binbackE[k,0,*]))
           y2 = double(transpose(binback[k,1,*]))
           y2err = double(transpose(binbackE[k,1,*]))
           yptitle='Back Ratio'
           if n_elements(timebin) NE 0 then begin
              print,"Binning not set up for individual background fluxes"
              return
           endif
           y = y / y2
        end
        keyword_set(background): begin
           y = double(transpose(binback[k,0,*]))
           yerr = double(transpose(binbackE[k,0,*]))
           y2 = double(transpose(binback[k,1,*]))
           y2err = double(transpose(binbackE[k,1,*]))
           yptitle='Flux (DN)'
           if n_elements(timebin) NE 0 then begin
              print,"Binning not set up for individual background fluxes"
              return
           endif
           individual=1
        end
        keyword_set(individual): begin
           y = double(transpose(binind[k,0,*]))
           yerr = double(transpose(binindE[k,0,*]))
           y2 = double(transpose(binind[k,1,*])) * reffactor
           y2err = double(transpose(binindE[k,1,*])) * reffactor
           yptitle='Flux (DN)'
           if n_elements(timebin) NE 0 then begin
              print,"Binning not set up for individual star fluxes"
              return
           endif
        end
        else: begin
           y = double(transpose(binfl[k,*]))
           yerr = double(transpose(binfle[k,*]))
           yptitle='Flux Ratio'
        end
     endcase

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
     ;; also for rejected state parameter points, remove those if
     ;; using a slit loss model
     goodp = where(finite(y) EQ 1 and badarray EQ 0)
     if goodp NE [-1] then begin
        y = y[goodp]
        yerr = yerr[goodp]
        if keyword_set(individual) then begin
           y2 = y2[goodp]
           y2err = y2err[goodp]
        endif
        tplot = tplot[goodp]
        airmass = airmass[goodp]
        thisUT = utgrid[goodp]
        if keyword_set(slitmod) then inputX = inputX[*,goodp]
        offp = where(tplot LT hstart OR tplot GT hend)
     endif else begin
        thisUT = utgrid
     endelse
     stdoff = stddev(y[offp])

     if keyword_set(slitmod) then begin
        smoothP = [0,1,2,3,5,6]
        nsmoothP = n_elements(smoothP)
        if n_elements(psmooth) EQ 0 then smoothsize = 30 else smoothsize=psmooth
        for i=0l,nsmoothP-1l do begin
           inputX[smoothP[i],*] = smooth(inputX[smoothP[i],*],smoothsize)
        endfor
     endif

     if not keyword_set(nonormalize) then begin
        yerr = yerr / median(y[offp])
        y = y / median(y[offp])
        if n_elements(y2) GT 0 then y2 = y2/median(y2[offp])
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
                 thisUT = thisUT[goodp]
                 if keyword_set(individual) then begin
                    y2 = y2[goodp]
                    y2err = y2err[goodp]
                 endif
                 airmass = airmass[goodp]
                 if keyword_set(slitmod) then inputX = inputX[*,goodp]
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
              thisUT = thisUT[goodp]
              if keyword_set(slitmod) then inputX = inputX[*,goodp]

              ;; fit result to a robust line
              rlinefit = robust_linefit(tplot,divbycurveclip1,yfit)
              ;; divide by the line to flatten out
              yflat = divbycurveclip1 / yfit
              rsigma = robust_sigma(yflat)
              maxval = max(yflat,maxp)
              minval = min(yflat,minp)
              nombysigma = (yflat - mean(yflat))/rsigma
              
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
                 thisUT = thisUT[goodp]
                 if keyword_set(slitmod) then inputX = inputX[*,goodp]
                 offp = where(tplot LT hstart OR tplot GT hend)
              endif
           endif
        endelse
        
        
     endif

     if n_elements(timebin) NE 0 then begin
        ;; Bin the series in time
        ybin = avg_series(tplot,y,y/yerr,timeGrid,tsizes,weighted=1,$
                  oreject=TsigRejCrit,eArr=yerrOut,/silent,errIn=yerr,stdevArr=stdevArr,$
                         median=tmedian)
        if n_elements(diffInd) NE 0 then begin
           if k EQ diffInd and keyword_set(differential) then begin
              ;; make sure the sigma rejection is off for the
              ;; normalization bin
              ybin = avg_series(tplot,y,y/yerr,timeGrid,tsizes,weighted=1,$
                                eArr=yerrOut,/silent,errIn=yerr,stdevArr=stdevArr,$
                               median=tmedian)
           endif
        endif

        airbin = avg_series(tplot,airmass,fltarr(n_elements(airmass)),timeGrid,tsizes,$
                            weighted=0,/silent,median=tmedian)
        utbin = avg_series(tplot,thisUT,fltarr(n_elements(thisUT)),timeGrid,tsizes,$
                           weighted=0,/silent,median=tmedian)
        tplot = tmiddle

        y = ybin
        fracPhotonarr[k] = median(yerrOut) / median(ybin)
;        yerr = yerrOut
        yerr = stdevArr
        airmass = airbin
        thisUT = utbin
        offp = where(tplot LT hstart OR tplot GT hend)

        ;; For any binned flux or error that are NAN, remove
        goodp = where(finite(y) EQ 1 and yerr GT 0)
        if goodp NE [-1] then begin
           y = y[goodp]
           yerr = yerr[goodp]
           if keyword_set(individual) then begin
              y2 = y2[goodp]
              y2err = y2err[goodp]
           endif
           tplot = tplot[goodp]
           airmass = airmass[goodp]
           thisUT = utgrid[goodp]
           if keyword_set(slitmod) then inputX = inputX[*,goodp]           
           offp = where(tplot LT hstart OR tplot GT hend)
        endif

        if keyword_set(offtranserr) then begin
           fitY2 = linfit(tplot[offp],y[offp])
           Offresid2 = y[offp] - (fitY2[0] + fitY2[1]*tplot[offp])
           if keyword_set(StraightupRMS) then begin
              rstdevOff2 = robust_sigma(y[offp])
           endif else rstdevOff2 = robust_sigma(Offresid2)
           yerr = fltarr(n_elements(yerr)) + rstdevOff2
        endif


     endif else fracPhotonArr[k] = median(yerr)/median(y)

     ;; Linearly detrend
     if keyword_set(lindetrend) then begin
        fitY = linfit(tplot[offp],y[offp])
        y = y/(fitY[0] + fitY[1]*tplot)
        if n_elements(y2) NE 0 then begin
           fitY2 = linfit(tplot[offp],y2[offp])
           y2 = y2/(fitY[0] + fitY[1]*tplot)
        endif
     endif

     if total(finite(y)) GT 2 and total(finite(yerr)) GT 2 then begin
        ;; if keyword set, replace the error w/ the off transit stddev

        ;find the range where 95% or more of the points are shown
        if keyword_set(individual) then begin
           ycomb = [y,y2] ;; combine both stars into one array
           sorty = sort(ycomb)
           ylength = n_elements(ycomb)
           ylowerL = ycomb[sorty[ceil(5E/100E*float(ylength))]] * 0.8
           yUpperL = ycomb[sorty[ceil(95E/100E*float(ylength))]] * 1.1
           ydynam = [ylowerL,yUpperL]
        endif else begin
           sorty = sort(y)
           ylength = n_elements(y)
           case 1 of
              keyword_set(fullrange): ydynam = [0,0]
;              keyword_set(oneprange): ydynam = [0.975,1.005]
              keyword_set(oneprange): ydynam = [0.99,1.005]
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

        if keyword_set(quadfit) then begin
           result = poly_fit(tplot[offp],y[offp],2,measure_errors=yerr[offp],yfit=yfit)
           Offresid = y[offp] - yfit
        endif else begin
           fitY = linfit(tplot[offp],y[offp])
           Offresid = y[offp] - (fitY[0] + fitY[1]*tplot[offp])
        endelse
        fracRMSarr[k] = robust_sigma(Offresid)/median(y[offp])

        ;; Show the off transit fit
;        oplot,tplot,fity[0] + fity[1]*tplot,color=mycol('red')
        if not keyword_set(noplots) then begin
           print,'Frac lin corr robust sigma for ',wavname[k],': ',fracRMSarr[k]

           if keyword_set(offtranserr) then begin
              
              if keyword_set(straightupRMS) then begin
                 rstdevOff = robust_sigma(y[offp])
              endif else rstdevOff = robust_sigma(Offresid)
              yerr = fltarr(n_elements(goodp)) + rstdevOff
           endif
           
           if keyword_set(psplot) and (not keyword_set(singleplot) OR k EQ 0)  then begin
              plotnmpre = 'plots/spec_t_series/tser_'+wavname[k]
              device,encapsulated=1, /helvetica,$
                     filename=plotnmpre+'.eps'
              if keyword_set(singleplot) then begin
                 device,xsize=PSSingleXsize, ysize=PSSingleYsize,decomposed=1,/color
              endif else device,xsize=PSplotXsize, ysize=PSplotYsize,decomposed=1,/color
           endif

           if keyword_set(singleplot) then begin
              if k EQ 0 then begin
                 ;; Set up everything for the first time plot
                 myNoerase=0
                 yptitle= yptitle + ' + Offset'
                 myXtitle='Orbital Phase'
                 tickformat=''
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
              if n_elements(custSep) EQ 0 then begin
                 if keyword_set(differential) then spacing = 0.08E else begin
                    spacing=separationA[0]
                 endelse
              endif else spacing = custSep
              ydynam=[1E - spacing * (1+Nwavbins),1+spacing]
              offset = k * spacing
              myTitle=''
           endif else begin
              myNoErase=0
              offset = 0
              myTitle=wavname[k]+' um Flux'
              tickformat=''
              myXtitle='Orbital Phase'
           endelse
           if keyword_set(custTitle) then begin
              myTitle = custTitle
           endif
           
           if keyword_set(custyrange) then ydynam = custyrange
           if (1 - keyword_set(singleplot)) OR k EQ 0 then begin
              ;; If in regular mode, plot each time, but in singleplot
              ;; mode, it only should draw axis the first time
              if keyword_set(showBJD) then begin
                 myXstyle = 8 + 1
;                 myPposition = [0.05,0.05,0.95,0.9]
                 moveTitle = myTitle
                 myTitle = ''
              endif else begin
                 myXstyle = 1
              endelse
              plot,tplot,y,psym=4,$
                   xtitle=myXtitle,$
                   title=myTitle,$
                   ytitle=yptitle,$
                   yrange=ydynam,ystyle=1,/nodata,xstyle=myXstyle,$
                   noerase=myNoErase,xthick=2.5,ythick=2.5,$
                   xtickformat=tickformat,ytickformat=tickformat,$
                   xrange=myXrange,xmargin=custxmargin,ymargin=custymargin
               if keyword_set(showBJD) then begin
                 refval = double(ev_round(min(utgrid),0.1D))
                 if min(utgrid - refval) LT -0.03 then refval = refval - 0.1D
;                 showUTrel = thisUT - refval
                 allphase = (utgrid - tmid)/ planetdat.period mod 1.0D
                 allphase = fold_phase(allphase,secondary=secondary)
                 UTlinfit = poly_fit(allphase,utgrid - refval,1)
;                 UTlinfit = poly_fit(tplot,thisUT - refval,1)
                 UTrelLimits = poly(!x.crange,transpose(UTlinfit))
                 axis,xaxis=1,xrange=UTrelLimits,xtitle='BJD - '+string(refval,format='(F9.1)'),$
                      xstyle=1
                 Xtitlepos = mean(!x.crange)
                 Ytitlepos = max(!y.crange) + 0.17E * (max(!y.crange) - min(!y.crange))
                 xyouts,Xtitlepos,Ytitlepos,moveTitle,alignment=0.5,charsize=0.7
                 print,'Range of utgrid = '
                 print,min(utgrid),format='(F16.5)'
                 print,max(utgrid),format='(F16.5)'
                 print,'Range used = '
                 print,UTrelLimits + refval,format='(F16.5)'
              endif
           endif
           if k mod 2 EQ 0 then dataColor=!p.color else dataColor=mycol('red')
           if not keyword_set(differential) then begin
              if keyword_set(showclipping) then begin
                 oplot,tplot,y,psym=5,color=mycol('red') ;; original data
                 oplot,tplotdivcurves,divbycurve,psym=4  ;; divided by light curve
                 oplot,tplot,yfit,color=mycol('blue')    ;; fitted line to curve
              endif else begin
                 if n_elements(timebin) EQ 0 then begin
                    oplot,tplot,y-offset,psym=4,color=dataColor
                 endif else begin
                    oploterror,tplot,y-offset,tsizes/2E,yerr,psym=3,$
                               hatlength=0,thick=2,color=dataColor
                 endelse
              endelse 
           endif
           if keyword_set(individual) then begin
              oplot,tplot,y2 - offset,psym=4,color=mycol('blue')
              al_legend,['Planet Host','Reference Star'],$
                        psym=[4,4],color=mycol(['black','blue']),$
                        /right,/clear
           endif
           if keyword_set(singleplot) and not keyword_set(skipwavl) then begin
              if valid_num(wavname[k]) then begin
                 wavelabel = wavname[k]+' '+cgGreek('mu')+'m'
              endif else begin
                 wavelabel = wavname[k]
              endelse
              if n_elements(wavelabelcsize) EQ 0 then begin
                 wavelabelcharS = 1.0
                 wavelabelcol = !p.color
              endif else begin
                 wavelabelcharS = wavelabelcsize
                 wavelabelcol = dataColor
              endelse
              wavelabelposx = !x.crange[1]-0.1*(!x.crange[1]-!x.crange[0])
              wavelabelposy = median(y)-offset
              xyouts,wavelabelposx,wavelabelposy,$
                     charsize=wavelabelcsize,$
                     [wavelabel],alignment=0.5,color=wavelabelcol
           endif
           ;; print the stdev for y for off points
;        print,'Fractional off transit Stdev in F for ',wavname,': ',stddev(y[offp])/mean(y[offp])
;        print,'Fractional off transit Robust sigma for ',wavname,': ',robust_sigma(y[offp])/median(y[offp])
           ;; try fitting the off transit to a function first
           
           ;; show the transit epochs
           drawy = [!y.crange[0],!y.crange[1]]
           oplot,[hstart,hstart],drawy,color=mycol('brown'),linestyle=2
           oplot,[hend,hend],drawy,color=mycol('brown'),linestyle=2
           
           ;; Show the jump point
           if keyword_set(showjump) then begin
              tjumpJD = date_conv('2013-08-15T09:39.00','J')
              tjump = (tjumpJD - tmid)/planetdat.period
              tjump = tjump mod 1D
              plots,[tjump,tjump],drawy,color=mycol('red'),thick=2,linestyle=2
           endif
           
           ;; show the off-transit fit for differential measurements (as
           if keyword_set(differential) then begin
              if not keyword_set(fitcurve) then begin
                 ;; long as fits aren't being performed)
                 
                 if keyword_set(quadfit) then begin
                    oplot,tplot,(result[0] + result[1] *tplot + result[2] * tplot^2)-offset,$
                          thick=2,color=mycol('blue')
                 endif else begin
                    oplot,tplot,(fitY[0] + fitY[1] *tplot)-offset,thick=2,color=mycol('blue')
                 endelse
              endif
              if n_elements(timebin) EQ 0 then tsizes = fltarr(n_elements(tplot)) + tplot[1]-tplot[0]
              oploterror,tplot,y-offset,tsizes/2E,yerr,psym=3,hatlength=0,thick=2,$
                         color=dataColor
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
                 plotnmpre = 'plots/error_distrib/error_hist_'+wavname[k]
                 device,encapsulated=1, /helvetica,$
                        filename=plotnmpre+'.eps'
                 device,xsize=PSplotXsize, ysize=PSplotYsize,decomposed=1,/color
                 
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
              if wavname[k] EQ 'z-prime' then begin
                 change_kernels,kernchoice[0]
              endif else begin
                 change_kernels,kernchoice[1]
              endelse
              nmcmcPars = n_elements(mcmcPars[0,*])
              mcmcShowP = 350
              phaseShow = findgen(mcmcshowP)/float(mcmcShowP) * (max(tplot)-min(tplot)) + min(tplot)
              ;; First find the mean function with expression evaluation
              meanfunctest = expression_eval(modelExpr,phaseShow,transpose(mcmcPars[k,0:8]))
                                ;oplot,phaseShow,meanfunctest-offset,color=mycol('blue'),thick=2
              ;; Find the residual vector
              meanfuncdat = expression_eval(modelExpr,tplot,transpose(mcmcPars[k,0:8]))
              mcmcResid = y - meanfuncdat
              ;; Find the inverse covariance matrix using the likelihood
              ;; function which also needs the inverse covariance matrix
              likeL = ev_leval([transpose(mcmcPars[k,9:nmcmcPars-1]),0],x=tplot,yin=y,yerr=yerr,Cinv=Cinv)
              mcmcModel = fltarr(n_elements(phaseShow))
              for m=0l,mcmcShowP-1l do begin
                 columnvec = cov_kernel(phaseShow[m] - tplot,transpose(mcmcpars[k,9:nmcmcPars-1]))
                 mcmcModel[m] = meanfunctest[m] + columnvec ## Cinv ## transpose(mcmcResid)
              endfor
              oplot,phaseShow,mcmcModel-offset,color=mycol('blue'),thick=2
              
              ;; Calculate a Chi-squared, even if it's not the
              ;; right thing
              mcmcModelmatch = fltarr(n_elements(y)) ;; match the # of points (1 model pt per data pt)
              for m=0l,n_elements(y)-1l do begin
                 columnvec = cov_kernel(tplot[m] - tplot,transpose(mcmcpars[k,9:nmcmcPars-1]))
                 mcmcModelmatch[m] = meanfuncdat[m] + columnvec ## Cinv ## transpose(mcmcResid)
              endfor
              pseudoresids = (mcmcModelmatch - y) ;; sort of like the residuals
              pseudochisquare = total((pseudoresids/yerr)^2)
              print,'Psuedo Chi-square = ',pseudochisquare
              if keyword_set(showNommcmc) then oplot,tplot,meanfuncdat-offset,color=mycol('green')
           endif
           
           if keyword_set(showKep) then begin
              ;; Show the Kepler Light Curve
              np = 1024E
              kphaseS = min(tplot) + findgen(round(np))/np * $
                        (max(tplot,/nan)-min(tplot,/nan))
              keplerF = kepler_func(kphaseS,1E)
              oplot,kphaseS,keplerF - offset,color=mycol('dgreen'),thick=4
              if keyword_set(labelKep) then begin
                 kscLabX = (!x.crange[1] - !x.crange[0]) * 0.1 + !x.crange[0]
                 kscLabY = (!y.crange[1] - !y.crange[0]) * 0.5 + !y.crange[0]
                 xyouts,kscLabX,kscLabY,'Kepler 13-Month Avg',color=mycol('dgreen'),$
                        charsize=0.6
              endif
           endif
           
           if keyword_set(fitcurve) then begin
              ;; fit the data curve
              start=double([planetdat.p,planetdat.b_impact,u1parm,u2parm,$
                            planetdat.a_o_rstar,1.0D,0D,0D,0D,0D,0D,0D])
              
              
              ;; Here's where you choose the model to fit to the data
              case 1 of
                 keyword_set(sinfit): begin
                    expr = 'P[0] * sin((X - P[2]) * 6.2831853D / P[1]) +  eval_legendre(X,P[5:10])'
                    start = fltarr(11)
                    start[0] = 0.002
                    start[1] = 0.175E
                    start[2] = -0.045E;-0.133E
                    start[3:4] = 0E
                    start[5] = 1E
                    start[6:10] = 0E
                 end
                 keyword_set(kepdiff): begin
                    expr = '(kepler_func(X,P[0]+1D) / kepler_func(X,1D)) *  eval_legendre(X,P[5:10])'
                 end
                 keyword_set(slitmod): begin
                    if keyword_set(individual) then begin
                       expr = 'vslit_approx(X[0,*] - P[5],P[8],X[2,*],X[5,*])'
                    endif else begin
                       expr = 'vslit_approx(X[0,*] - P[5],P[8],X[2,*],X[5,*])/'+$
                              'vslit_approx(X[1,*] - P[5] - P[6],P[8],X[3,*] * P[7],X[6,*])'
                    endelse
                    start[5] = 0.0E ;; start the position as 0
                    start[7] = 1.0E ;; start the relative FWHM ratio as 1.0
                    if strmatch(usedate,'*2012*') OR strmatch(usedate,'*2013*') OR $
                       strmatch(usedate,'*2011*') then begin
                       H = 20.49/2E     ;; for the Aladin Detector
                    endif else H = 30.5/2E ;; for the Hawaii II detector
                    start[8] = H
                    start[10] = 1E  ;; start the offset parameter as 1.0
                    start = [start,0E] ;; last Legendre parameter should be 0E as a start
                    case 1 of
                       keyword_set(secondary): begin
                          start[0] = 0.002E
                          expr = expr+' * sec_eclipse(X[4,*]-P[9],P[0],P[1],P[4]) * eval_legendre(X,P[10:12])'
                       end
                       keyword_set(kepfit): begin
                          expr = expr+' * kepler_func(X[4,*],P[0]) * eval_legendre(X,P[10:12])'
                          start[0] = 1E
                       end
                       else: begin
                          expr = expr+' * quadlc(X[4,*]-P[9],P[0],P[1],P[2],P[3],P[4])* eval_legendre(X,P[10:12])'
                       end
                    endcase
                 end
                 keyword_set(differential): begin
                    if keyword_set(secondary) then begin
                       expr = 'sec_eclipse(X-P[11],P[0],P[1],P[4]) * eval_legendre(X,P[5:10])'
                       start[0] = 0.002E
                    endif else begin
                       quadlcArg ='X'
                       for m=0l,4 do begin
                          quadlcArg=quadlcArg+','+strtrim(start[m],1)
                       endfor
                       expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* (P[5] + X * P[6] + X^2 * P[7] + X^3'+$
                              ' * P[8])/quadlc('+quadlcArg+')'
                    endelse
                    
                 end
                 keyword_set(kepfit): begin
                    expr = 'kepler_func(X-P[11],P[0]) *  eval_legendre(X,P[5:10])'
;                 expr = 'parameterized_kep(X-P[11],P[0]) *  eval_legendre(X,P[5:10])'
                 end
                 keyword_set(secondary): begin
                    expr = 'sec_eclipse(X-P[11],P[0],P[1],P[4]) * eval_legendre(X,P[5:10])'
                    start[0] = 0.002E
                 end
                 else: begin
                    expr = 'quadlc(X-P[11],P[0],P[1],P[2],P[3],P[4])* eval_legendre(X,P[5:10])'
                 end
              endcase
              
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Here's the long parameter setting section - fixed, free or limited
              
              ;; By default all parameters start out as fixed and >=0
              pi = replicate({fixed:1, limited:[1,0], limits:[0.0E,0.0E]},nparams)
              
              ;; make sure the Rp/R* parameter is free
              if not keyword_set(fixrad) then pi[0].fixed = 0 
              
              ;; Limb Darkening parameter adjustment (3 options)
              if keyword_set(freelimbquad) then begin
                 ;; if asked to, free the quadratic limb darkening parameter
                 pi[2].fixed = 0
                 pi[3].fixed = 0
              endif
              if keyword_set(freelimblin) then pi[2].fixed = 0
              ;; Let the limb darkening be + or -
              pi[2].limited = [0,0]
              pi[3].limited = [0,0]
              
              ;; Here are a few modifications to the Rp/R* parameter
              ;; (letting it be <0 for differential Rp/R* fitting)
              if not keyword_set(kepfit) and $
                 not keyword_set(kepdiff) and $
                 not keyword_set(secondary) then begin
                 pi[0].limited = [1,1] ;; make sure Rp/R* is limited
                 pi[0].limits = [0D,1D] ;; Keep Rp/R* between 0 and 1
              endif else begin
                 pi[0].limited = [0,0]
                 pi[0].limits = [0D,0D]
              endelse


              ;; Transit epoch fitting
              if keyword_set(fitepoch) then begin
                 pi[11].fixed = 0
                 pi[11].limited = [1,1]
;              pi[11].limits = [-0.05,0.05]
                 pi[11].limits = [-0.1,0.1]
                 start[11] = 0.00
              endif
              
              ;; Details of slit model adjustments (making them free)
              if keyword_set(slitmod) then begin
                 pi[5:7].fixed=0      ;; free the slit model params
                 pi[5:6].limited = [1,1]            ;; limit the slit position from going to far
                 pi[5:6].limits = [-start[8],start[8]] ;; slit half width
                 if keyword_set(individual) then begin
                    pi[6:7].fixed=1 ;; only look at 1 star at a time
                 endif
              endif
              
              if keyword_set(sinfit) then begin
                 pi[0:2].fixed =0  ;; allow amplitude, phase and period be free
                 pi[3:4].fixed =1 ;; make sure unused parameters are fixed
                 ;;pi[5:10] will be handled below by the legendre
                 ;;stuff
                 pi[0].limited = [1,0]
;                 pi[0].limited = [1,0]
                 pi[0].limits = [0,0] ;; ensure amplitude is positive
                 pi[1].fixed = 1
;                 pi[1].limited = [1,1]
;                 pi[1].limits = [0.16,0.19]
;                 pi[1].limits = [1E-4,0] ;; ensure you don't divide by zero for period
                 ;pi[2].limited = [1,1]
                 pi[2].limits = [-0.15,-0.12]
                 ;pi[2].limits = [0,0.2]
              endif

              ;; Here's where the polynomial terms are adjusted
              if keyword_set(slitmod) then begin
                 startLeg = 10 ;; first legendre polynomial coefficient
                 lastLeg = 12  ;; last Legendre polynomial in model
              endif else begin
                 startLeg = 5 ;; first legendre polynomial coefficient
                 lastLeg = 10 ;; last Legendre polynomial in model
              endelse
              if n_elements(legOrder) EQ 0 then legOrder =1 ;; the default is a linear baseline
              assert,lastLeg,'>=',legOrder+startLeg,"More polynomial terms are asked for than allowed by model"
              if keyword_set(fixoffset) then begin
                 ;; nothing to do b/c the rest of the parameters are already fixed
              endif else begin
                 endLeg = startLeg + LegOrder ;; final Legendre polynomial
                 pi[startLeg:endLeg].fixed = 0   ;; free them!
                 pi[startLeg:endLeg].limited = [0,0] ;; no limits
              endelse
              
              ;; The freeall and fixall keyword override the rest of those
              if keyword_set(freeall) then begin
                 pi[*].fixed = 0 ;; free all parameters if asked to
              endif 
              if keyword_set(fixall) then begin
                 pi[*].fixed = 1 ;; fix all parameters if asked to
              endif
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
              ;; Here's where a Levenberg-Marquardt fit is applite to the data
              if keyword_set(slitmod) then begin
                 result = mpfitexpr(expr,inputX,y,yerr,start,parinfo=pi,perr=punct,/quiet)
                 modelY = expression_eval(expr,inputX,result)
                 modelY1 = modelY
                 modelX = tplot
              endif else begin
                 if keyword_set(boot) then begin
                    result = boot_mp(expr,tplot,y,yerr,start,parinfo=pi,perr=punct,/quiet)
                 endif else begin
                    result = mpfitexpr(expr,tplot,y,yerr,start,parinfo=pi,perr=punct,/quiet)
                 endelse
                 modelPts = 512l
                 ntpoints = n_elements(tplot)
                 modelY = expression_eval(expr,tplot,result)
                 modelX = findgen(modelPts)/(modelPts-1l) * (tplot[ntpoints-1] - min(tplot,/nan)) + min(tplot,/nan)
                 modelY1 = expression_eval(expr,modelX,result)
              endelse
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
                 oplot,modelX,modelY2 - offset,color=mycol('orange'),thick=2
                 al_legend,['Best-Fit Radius','Nominal Radius'],$
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
                    plotnmpre = 'plots/residual_series/residuals_'+wavname[k]
                    device,encapsulated=1, /helvetica,$
                           filename=plotnmpre+'.eps'
                    device,xsize=PSplotXsize, ysize=PSplotYsize,decomposed=1,/color
                    
                 endif
                 
                 ylowerL = resid[sorty[ceil(5E/100E*float(ylength))]]
                 yUpperL = resid[sorty[floor(95E/100E*float(ylength))]]
                 ydynam = [-1E,1E] * max(abs([ylowerL,yUpperL])) * 4E
                 
                 if n_elements(custresidYrange) NE 0 then ydynam = custresidYrange
                 
                 overplotMarg = [13,14]
                 plot,tplot,resid,yrange=ydynam,$
                      title='Residuals at '+wavname[k],$
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
                 print,'RMS Residuals (%) for '+wavname[k],'um',(stddev(y - modelY))/median(y)*100E
              endif
           endif
           if not keyword_set(fitcurve) then begin
              ntplot = n_elements(tplot)
              modelY = fltarr(ntplot)
              resid = fltarr(ntplot)
           endif
           if keyword_set(showmcmc) then begin
              if keyword_set(useGPasfit) then begin
                 resid = pseudoresids * 100E
              endif else resid = mcmcResid * 100E
              modelY = meanfuncdat
           endif
           
           cleanTSerName = 'data/cleaned_tim_ser/timeser_'+wavname[k]+'um_.txt'
           forprint,tplot,y,yerr,modelY,resid,$
                    textout=cleanTSerName,$
                    comment='#Phase  Flux  Fl_err  Model_fl   Residual for '+wavname[k]+'um (%)',$
                    /silent
           if n_elements(cleanlist) EQ 0 then begin
              cleanlist = cleanTSerName
           endif else cleanlist = [cleanlist,cleanTSerName]
           
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
     endif
     
  endfor
  
  ;; save the RMS of off transit fluxu
  ;; If there's no wavelength binning
  if n_elements(timebin) EQ 0 then tsizes = [0]
  save,bingrid,fracRMSarr,tsizes,bingridmiddle,binsizes,$
       planetdat,fracPhotonarr,$
       filename='data/rmsdata.sav'
  
  ;; Save the clean list name
  if n_elements(cleanlist) NE 0 then begin
     forprint,/nocomment,cleanlist,textout='data/cleaned_list.txt',format='(A)',/silent
  endif
  
  if keyword_set(fitcurve) then begin
     ev_print_params,wavname,paramnames,resultarr,resultarrE,pi,k,/skipsig
     
     ;; save the radius data
     forprint,bingridmiddle[*],binsizes,plrad,plrade,$
              textout='radius_vs_wavelength/radius_vs_wavl.txt',$
              comment='#Wavelength(um) Binsize (um)  Rp/R*   Rp/R* Error',/silent
     
     depth = plrad^2 * 1E6
     depthErrP = (plrad + plrade)^2 * 1E6 - depth
     depthErrM = depth - (plrad - plrade)^2 * 1E6 
     
     forprint,bingridmiddle[*],binsizes,depth,depthErrP,depthErrM,$
              textout='radius_vs_wavelength/depth_vs_wavl.txt',$
              comment='#Wavelength(um) Binsize (um)  Depth (ppm)   Depth Error+(ppm), Depth Error-(ppm)',/silent
     
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
  endif
  
  if keyword_set(singleplot) and not keyword_set(skipReset) then begin
     !p.multi = 0
     !Y.Omargin = [0,0]
  endif
  
  if keyword_set(psplot) then begin
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
end
