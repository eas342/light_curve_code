pro simulate_series
;; Makes a random time series with correlated noise

  Npoints=600l
;  Npoints=1000l
  X = findgen(Npoints)
;  RandomSer = randomn(0,Npoints)

  Nrealizations = 5l
  
  startseeds = [0,2,3,5,7]
  colarray = mycol(['red','green','blue','orange','yellow'])
  for j=0l,Nrealizations-1l do begin 

     RandomSer = randomn(startseeds[j],Npoints)

     Y = fltarr(Npoints) ;; blank array to start

     Nlength = 99l
     Nstrength = 0.013E
;     Nlength = 20l
;     Nstrength = 0.0505E
  ;; Initialize the first Nlength points
     Y[0l:Nlength-1l] = RandomSer[0l:Nlength-1l]
     MixingA = fltarr(Nlength) + Nstrength
     for i=Nlength,Npoints-1l do begin
        Y[i] = RandomSer[i] + total(MixingA[0l:Nlength-1l] * Y[i-Nlength:i-1])

;     Y[i] = RandomSer[i]*(1E - Mixing1 - Mixing2 - Mixing3) + Y[i-1] * Mixing1 + $
;            Y[i-2] * Mixing2 + Y[i-3] * Mixing3
     endfor
     if j EQ 0l then begin
        plot,x,y,yrange=[-10,10]
     endif else begin
        oplot,x,y,linestyle=j,color=colarray[j]
     endelse

  endfor

  
  steparray = lindgen(Npoints)
  autoC = a_correlate(y,steparray)

;  plot,steparray[1l:Npoints-1l],autoC[1l:Npoints-1l],/xlog
;  stop

end
