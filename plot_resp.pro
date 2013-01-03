pro plot_resp,psplot=psplot,noskip=noskip,custxrange=custxrange,$
              custyrange=custyrange,individual=individual
;; the psplot keyword tells the program to make a Postscript file
;; noskip does not skip the first few files in the series (detector
;; ramp up)
;; custxrange allows a custom x range instead of the default
;; custyrange allows a custom y range instead of the default
;; individual -- allows you to look at the stars one at a time

  restore,'data/fluxsave.sav'
  ;; should pick up flarr(average pixel flux per file and box),
  ;;t(time),expt(seconds),nbox(number of
  ;; extraction boxes) and nfile(number of files),flerr(error in
  ;; average pixel flux),wavl(the approximate wavelength for each box)
  ;; t (time in JD for each frame)

  ;; convert expt from a string to a float
  expt = float(expt)

  ;; look up which date the data come from
  JDates = [date_conv('2011-12-23T08:06:00.00','JD'),$
            date_conv('2011-12-29T08:06:00.00','JD'),$
            date_conv('2012-01-04T08:06:00.00','JD')]
  nightnames = ['dec_23','dec_29','jan_04']
  tabinv,JDates,t[0],night
  night = round(night)

  ;; get the transit times
  readcol,'transit_info/'+nightnames[night]+'_t_time.txt',skipline=1,$
          format='(A,A)',epoch,tepoch
  t1stcontact = date_conv(tepoch[0],'JULIAN')
  t4thcontact = date_conv(tepoch[1],'JULIAN')
  tmidtrans = (t1stcontact + t4thcontact)/2D

  trel = (t - tmidtrans)*(24E*60E) ;; time relative to mid-transit in minutes

  if keyword_set(psplot) then begin
     if keyword_set(noskip) then begin
        plotnm = 'plots/basic_transit_noskip.eps'
     endif else begin
        plotnm = 'plots/basic_transit.eps'
     endelse
     set_plot,'ps'
     !p.font=0
     device,encapsulated=1, /helvetica,$
            filename=plotnm
     device,xsize=14,ysize=10,decomposed=1,/color
  endif

;  if keyword_set(noskip) then begin
  nfile = n_elements(flarr)
  set1 = lindgen(nfile)
;     set1 = where(expt EQ 10.0)
;  endif else begin
     ;; skip the first few frames
;     set1 = where(expt EQ 10.0 and lindgen(nfile) NE lindgen(6))
;  endelse
  x = trel[set1]
  if keyword_set(individual) then begin
     y = flarr[*,0] / mean(flarr[*,0]) ;; show the stars individually
     yerr = flerr[*,0] / mean(flarr[*,0])
     y2 = flarr[*,1] / mean(flarr[*,1])
     y2err = flerr[*,1] / mean(flarr[*,1])
  endif else begin
     y = flarr[*,0] / flarr[*,1]
     yerr = sqrt((flerr[*,0]/flarr[*,0])^2 + (flerr[*,1]/flarr[*,1])^2 )*y
  endelse

  ;; for now, don't show any error for time
  xerr = fltarr(n_elements(set1))
  
  if n_elements(custxrange) EQ 0 then custxrange = [0,0] ;; use default range (if not set)
  if n_elements(custyrange) EQ 0 then begin
     custystyle = 16
     custyrange = [0,0] ;; use default range (if not set)
  endif else custystyle = 1

  cgplot,x,y,psym=2,$
         xtitle='Time - Transit (min)',$
         ytitle='Flux Ratio',$
         ystyle=custystyle,/nodata,xrange=custxrange,$
         yrange=custyrange
  oploterror,x,y,xerr,yerr,psym=2

  if keyword_set(individual) then begin
     oploterror,x,y2,xerr,y2err,color=mycol('blue')
  endif

  ;; show the transit times
  plots,[t4thcontact - tmidtrans,t4thcontact - tmidtrans]*60E*24E,$
        !y.crange,color=mycol('blue')
  plots,[t1stcontact - tmidtrans,t1stcontact - tmidtrans]*60E*24E,$
        !y.crange,color=mycol('blue')


  print,'Stdev of Frame counts = ',stddev(y)/mean(y)*100E,'%'
  print,'Photon error of frame counts = ',mean(yerr)/mean(y)*100E,'%'
  if keyword_set(psplot) then begin
     device,/close
     cgPS2PDF,plotnm
     ;; return the plotting to normal
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
end
