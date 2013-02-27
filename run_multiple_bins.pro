pro run_multiple_bins
  night=['dec23','jan04']
  for i=0,n_elements(night)-1 do begin
     for j=5,5 do begin ;; cycle through # of wavelength bins
        case night[i] of
           'dec23': compile_spec,/optimal,nwavbins=j,/dec23
           'jan04': compile_spec,/optimal,nwavbins=j
        endcase
        plot_tim_ser,/fitrad,/freelimb,/quadfit,/offtranserr
        plot_rad_vs_wavl,/psplot
        spawn,'mv plots/rad_vs_wavl.png plots/rad_vs_wavl_'+night[i]+'_'+strtrim(j,1)+'pts.png'
        spawn,'mv plots/rad_vs_wavl.pdf plots/rad_vs_wavl_'+night[i]+'_'+strtrim(j,1)+'pts.pdf'
     endfor
  endfor
end
