pro run_multiple_bins
  for i=5,10 do begin
     compile_spec,/optimal,nwavbins=i
     plot_tim_ser,/fitrad,/freelimb,/quadfit
     plot_rad_vs_wavl,/psplot
     spawn,'mv plots/rad_vs_wavl.png plots/rad_vs_wavl_'+strtrim(i,1)+'pts.png'
  endfor
end
