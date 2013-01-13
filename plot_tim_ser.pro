pro plot_tim_ser,fitrad=fitrad,fitpoly=fitpoly,usepoly=usepoly,npoly=npoly,$
                 fullrange=fullrange,smartbin=smartbin,oneprange=oneprange,$
                 offtranserr=offtranserr,freelimb=freelimb,clarlimb=clarlimb
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


  ;; set the plot
  set_plot,'ps'
  !p.font=0

  ;; get the binned data
  restore,'data/specdata.sav'

  bingridmiddle = bingrid + binsizes/2E

  nut = n_elements(utgrid)
  nbin = Nwabins

  ;; get the transit times
  readcol,'transit_info/jan_04_t_time.txt',epoch,tepoch,format='(A,A)',$
          skipline=1

  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,F)',$
          skipline=1
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor

  u1parm = 0.0774E              ;HD 189733
  u2parm = 0.3097E

  tstart = date_conv(tepoch[0],'JULIAN')
  tend = date_conv(tepoch[1],'JULIAN')
  tmid = (tend + tstart)/2D

  ;; radius of a planet as a function of wavelength
  plrad = fltarr(nbin)*!values.f_nan
  plrade = fltarr(nbin)*!values.f_nan

  ;; orbital phase
  tplot = (utgrid - tmid)/planetdat.period
     
  ;; calculate start and end
  hstart = (tstart - tmid)/planetdat.period
  hend = (tend - tmid)/planetdat.period

  offp = where(tplot LT hstart OR tplot GT hend)

  ;; plot the binned data for aperture 2

  for k=0l,nbin-1l do begin
     specpt = k
     ;; choose the data from the normalization order
     if keyword_set(smartbin) then begin
        y = sbinfl[specpt,Nord,*,apPlot]
        yerr = sbinfle[specpt,Nord,*,apPlot]
     endif else begin
        y = binfl[specpt,Nord,*,apPlot]
        yerr = binfle[specpt,Nord,*,apPlot]
     endelse

     
     wavname = string(bingridmiddle[k,Nord],format='(F4.2)')

     if total(finite(y)) GT 0 and total(yerr) GT 0.0 then begin

        ;; if asked to, use the polynomial to correct the light curve
        ;; (only do it for valid polynomial fits)
        if keyword_set(usepoly) then begin
           if finite(ppolystruct.(0)[1,k]) then begin
              baseline = fltarr(nut)
              for i=0l,npoly-1l do begin
                 baseline = baseline + ppolystruct.(0)[i+1,k+1] * tplot^i
              endfor
              y = y - baseline
           endif
        endif

        ;; if keyword set, replace the error w/ the off transit stddev
        if keyword_set(offtranserr) then begin
           stdevOff = stddev(y[offp])
           yerr = fltarr(nut) + stdevOff
        endif

        ;find the range where 95% or more of the plots are shown
        sorty = sort(y)
        ylength = n_elements(y)
        case 1 of
           keyword_set(fullrange): ydynam = [0,0]
           keyword_set(oneprange): ydynam = [-0.01,0.01]
           else: ydynam = $
              y[sorty[[ceil(5E/100E*float(ylength)),floor(95E/100E*float(ylength))]]]
        endcase
        
        device,encapsulated=1, /helvetica,$
               filename='plots/spec_t_series/tser_'+wavname+'.eps'
        device,xsize=14, ysize=10,decomposed=1,/color
        plot,tplot,y,psym=2,$
             xtitle='Orbital Phase',$
             title=wavname+' um Flux / '+normtext+' um Flux',$
             ytitle='Flux Ratio - Mean',$
             yrange=ydynam
        oploterr,tplot,y,yerr

        ;; print the stdev for y
        print,'Stdev in F for ',wavname,': ',stddev(y)

        ;; show the transit epochs
        drawy = [-0.002,0.002]
        plots,[hstart,hstart],drawy,color=mycol('brown'),linestyle=2
        plots,[hend,hend],drawy,color=mycol('brown'),linestyle=2

        if keyword_set(fitrad) then begin
           ;; fit the data curve
           expr1 = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])'
           expr2 = ' - quadlc(X,'+strtrim(planetdat.p,1)+$
                   ','+strtrim(planetdat.b_impact,1)+$
                   ','+strtrim(u1parm,1)+$
                   ','+strtrim(u2parm,1)+$
                   ','+strtrim(planetdat.a_o_rstar,1)+') + P[5]'
           expr = expr1 + expr2
           pi = replicate({fixed:1, limited:[1,0], limits:[0.0E,0.0E]},6)
           pi[0].fixed = 0 ;; make sure the Rp/R* parameter is free
           pi[5].fixed = 0 ;; make sure the flux ratio offset is free
           ;; fix the impact parameter, limb darkening and AoR*
           if keyword_set(freelimb) then begin
              ;; if asked to, free the limb darkening parameter
              pi[2].fixed = 0
              pi[3].fixed = 0
           endif
           if keyword_set(clarlimb) then begin
              start=[planetdat.p,planetdat.b_impact,u1bin[k],u2bin[k],planetdat.a_o_rstar,0.0E]
           endif else begin
              start=[planetdat.p,planetdat.b_impact,u1parm,u2parm,planetdat.a_o_rstar,0.0E]
           endelse
           result = mpfitexpr(expr,tplot,y,yerr,start,parinfo=pi,perr=punct)
           
           ;; show the fit
           showphase = 0.08E
           nshowpts = 1024      ; number of points to show
           phtest = dindgen(nshowpts)/float(nshowpts)*showphase - showphase/2E
           ytest1 = quadlc(phtest,result[0],result[1],result[2],result[3],result[4],planetdat.a_o_rstar)
           ytest2 = quadlc(phtest,planetdat.p,planetdat.b_impact,u1parm,u2parm,planetdat.a_o_rstar)
           ytest = ytest1 - ytest2 + result[5]
           oplot,phtest,ytest,color=mycol('orange')

           ;; save the planet radius
           plrad[k] = result[0]
           plrade[k] = punct[0]
           

        endif
        if keyword_set(fitpoly) then begin
           ;; fit the data curve to a polynomial
           result = poly_fit(tplot,y[0l:nut-1l],npoly-1,measure_errors=yerr[0l:nut-1l],yfit=yfit)

           ;; show the fit
           showphase = 0.08E
           nshowpts = 1024      ; number of points to show
           ytest = fltarr(nshowpts)
           phtest = dindgen(nshowpts)/float(nshowpts)*showphase - showphase/2E
           for i=0,npoly-1 do begin
              ytest = ytest + phtest^i * result[i]
           endfor
           oplot,phtest,ytest,color=mycol('orange')
           
           polyparms[k,*] = result
        endif

        device, /close
        device,decomposed=0
        cgPS2PDF,'plots/spec_t_series/tser_'+wavname+'.eps',$
                  /delete_ps
           
        
     endif


  endfor
  
  ;; save the radius data
  
  forprint,bingridmiddle[*,Nord],plrad,plrade,$
           textout='radius_vs_wavelength/radius_vs_wavl_ord_'+ordernames[Nord]+'.txt',$
           comment='#Wavelength  Rp/R*   Rp/R* Error'

  if keyword_set(fitpoly) then begin
     ;; save the polynomial fit data
     openw,1,'radius_vs_wavelength/polynomial_params_ord_'+ordernames[Nord]+'.txt',$
           width=200
     printf,1,'#Wavelength ',format='($,A)'
     for i=0,npoly-1 do begin
        printf,1,'P[',strtrim(i,1),']   ',format='($,A,A,A)'
     endfor
     printf,1,'' ; new line
     for k=0,nbin-1l do begin
        printf,1,bingridmiddle[k,Nord],' ',format='($,A,A)'
        for i=0,npoly-1l do begin
           printf,1,polyparms[k,i],' ',format='($,A,A)'
        endfor
        printf,1,''
     endfor
     
;     forprint,bingrid[*,Nord],polyparms[*,0],polyparms[*,1],polyparms[*,2],$
;              polyparms[*,3],polyparms[*,4],$
;              textout='radius_vs_wavelength/polynomial_params_ord_'+ordernames[Nord]+'.txt',$
;              comment='#Wavelength P[0]      P[1]      P[2]      P[3]    P[4]'
     close,1
  endif
  
  device,decomposed=0
  set_plot,'x'
  !p.font=-1

end
