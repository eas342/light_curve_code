pro simulate_series
;; Makes a random time series with correlated noise

  Npoints=600l
  X = findgen(Npoints)
  RandomSer = randomn(0,Npoints)
  Y = fltarr(Npoints) ;; blank array to start
  Y[0] = RandomSer[0]
  Y[1] = RandomSer[1]
  Mixing1 = 0.8
  Mixing2 = 0.1
  for i=2l,Npoints-1l do begin
     Y[i] = RandomSer[i]*(1E - Mixing1 - Mixing2) + Y[i-1] * Mixing1 + Y[i-2] * Mixing2
  endfor

;  plot,x,y
;  wait,1
  
  steparray = lindgen(Npoints)
  autoC = a_correlate(y,steparray)

  plot,steparray[1l:Npoints-1l],autoC[1l:Npoints-1l],/xlog
;  stop

end
