pro print_times
  ;; Prints the transit times

;;  T0 = 2454159.4532D OLD Barge 2008 paper
  T0 = 2454159.452879D  ;; Bean paper 2009
  Tobs = date_conv('2012-01-04T10:00:00.00','JULIAN')
;;;  PeriodP = 1.5089557D Old Barge 2008 paper
  PeriodP = 1.5089656D 
  np = round( (Tobs - T0)/PeriodP)
  Tmid = T0 + double(np) * Periodp
;  PeriodUnc = 0.0000064D
  PeriodUnc = 0.000006D
  Tunc = double(np) * PeriodUnc
  print,'TMID (BJD) = ',Tmid,' +/- ',Tunc,format='(A,F,A,F)'

  ;; Change from BJD to UTC
  ;; get the planet coordinates
  readcol,'../param_input/planet_host_coordinates.txt',descrip,coord,format='(A,F)',$
          skipline=1
  RA = coord[0]
  Dec = coord[1]
;  print,bjd2utc(Tmid,RA,DEC),' second delay'
  print,'TMID = ',date_conv(Tmid,'FITS')

  duration = 0.104D ;; from exoplanets.org
  Tstart = Tmid - duration/2D
  Tend = Tmid + duration/2D
  print,'Tstart = ',date_conv(Tstart,'FITS')
  print,'Tend = ',date_conv(Tend,'FITS')


end
