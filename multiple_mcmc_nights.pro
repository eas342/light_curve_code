pro multiple_mcmc_nights

  ;; Make light curves, and do MCMC fitting for multiple nights
  ;; Make sure to empty the Cleaned_tim_ser folder first

  read,'Have you emptied the cleaned_time_ser folder?',junk

  read,'Have you emptied the param_unc folder?',junk

  for i=0l,1l do begin
     if i EQ 0 then begin
        compile_spec,/dec23,nwavbins=4
     endif else compile_spec,nwavbins=4
     plot_tim_ser,timebin=50
     try_mcmc
     gather_mcmc_radii
     if i EQ 0 then begin
        spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
              'radius_vs_wavelength/straightdec23_leg01fit_freelimblin_mcmc_hypers_'+$
              'free_free_yerr_004pt_absexp_kern.txt'
     endif else begin
        spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
              'radius_vs_wavelength/straightjan04_leg01fit_freelimblin_mcmc_hypers_'+$
              'free_free_yerr_004pt_absexp_kern.txt'
     endelse
  endfor


end
