pro remove_linear,utgrid,fluxArray,lamgrid
;; Takes a data array array with dimensions (#Wavelengths,Data), a
;; time array and wavelength array and de-trends all data along the
;; time direction
  
  ;;throw away all n_sigma events before de-trending
  firstCutSig = 12E
  
  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
          skipline=1
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor
  
  ;; get the transit times
  readcol,'transit_info/transit_epoch.txt',epoch,tepoch,format='(A,D)',$
          skipline=1
  tstart = tepoch[0]
  tend = tepoch[1]
  tmid = tepoch[2]
  
  ;; orbital phase
  tplot = fold_phase((utgrid - tmid)/planetdat.period)
  
  ;; calculate start and end
  hstart = fold_phase((tstart - tmid)/planetdat.period)
  hend = fold_phase((tend - tmid)/planetdat.period)
  
;     restore,'data/timedata.sav'
  nwavs = n_elements(lamgrid)
  offp = where(tplot LT hstart OR tplot GT hend)
  for i=0l,nwavs-1l do begin
     if total(finite(fluxArray[i,0,offp])) EQ 0 then goodp = [-1] else begin
        rstdoff = robust_sigma(fluxArray[i,0,offp])
        medoff = median(fluxArray[i,offp])
        
        goodp = where(abs(fluxArray[i,0,*] - medoff) LE firstCutSig * rstdoff and $
                      (tplot LT hstart OR tplot GT hend),complement=throwaways)
     endelse
     if n_elements(goodp) GT 5 then begin
        ytemp = fluxArray[i,0,goodp]
        
        yerr = fluxArray[i,0,goodp]
        tplottemp = tplot[goodp]
        
        ;; fit result to a robust line
        rlinefit = robust_linefit(tplottemp,ytemp)
        fluxArray[i,0,*] = fluxArray[i,0,*] / (rlinefit[0] + tplot * rlinefit[1])
     endif
  endfor

end
