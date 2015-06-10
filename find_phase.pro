function find_phase,utgrid
;; Returns the orbital phase (if given UTC_JDB)

  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
          skipline=1
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor

  nwavs = n_elements(lamgrid)
  ntimep = n_elements(utgrid)

  ;; Check if fluxarray has 3 dimensions or 2
  sz = size(fluxArray)
  ndim = sz[0]
  
  ;; get the transit times
  readcol,'transit_info/transit_epoch.txt',epoch,tepoch,format='(A,D)',$
          skipline=1
  tstart = tepoch[0]
  tend = tepoch[1]
  tmid = tepoch[2]
  
  ;; orbital phase
  tplot = fold_phase((utgrid - tmid)/planetdat.period)

  return,tplot

end
