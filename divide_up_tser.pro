pro divide_up_tser,initialize=initialize
;; Divides up the time series wavelengths into groups

  wholeRange = [0.85,2.44]
  wholeLength = wholeRange[1] - wholeRange[0]
  ngroups= 3
  nSubWaves = 9
  stepSize = wholeLength / ngroups
  nAllWavs = nSubWaves * ngroups
  gstarts = findgen(ngroups) * stepSize + wholeRange[0]
  gends = gstarts + stepSize
  
  for i=0l,ngroups-1l do begin
     thisWaveRange = [gstarts[i],gends[i]]
     compile_spec,/readc,/specsh,nwavbins=nSubWaves,custrange=thisWaveRange
     plot_tim_ser,/singlep,/fitcurve,/sinfit,/offtranserr,/lind,/psplot
  endfor
  compile_spec,/readc,/specsh,nwavbins=nAllWavs,custrange=wholeRange
  plot_tim_ser,/singlep,/fitcurve,/sinfit,/offtranserr,/lind

end
