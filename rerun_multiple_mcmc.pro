pro rerun_multiple_mcmc,docov=docov,both=both
;; Takes the MCMC chain data and re-runs portions of it
  ;; Make light curves, and do MCMC fitting for multiple nights
  ;; Make sure to empty the Cleaned_tim_ser folder first

;; docov -- makes covariance plots for the parameters and saves the
;;          eps files
;; both -- do both photometry and spectroscopy

  read,'Have you emptied the cleaned_time_ser folder?',junk

  read,'Have you emptied the param_unc folder?',junk

  discardpoints = 1000

  for i=0l,2l do begin
     if keyword_set(both) then begin
        case i of
           0: compile_both,/dec23
           1: compile_both ;jan 04
           2: compile_both,/dec29
        endcase
        radcode = 'both'
     endif else begin
        case i of
           0: compile_spec,/dec23,nwavbins=9
           1: compile_spec,nwavbins=9 ;; Jan 04
           2: compile_spec,/dec29,nwavbins=9
        endcase
        radcode = 'speX'
     endelse
     if i LE 1 then Npoints = 100 else Npoints = 50

     plot_tim_ser,timebin=Npoints,/offtranserr
     case i of
        0: datename = 'dec23'
        1: datename = 'jan04'
        2: datename = 'dec29'
     endcase
     cd,c=currentd
     fileopt = file_search(currentd+'/data/cleaned_tim_ser/*.txt')
     totfiles = n_elements(fileopt)
     for j=0l,n_elements(fileopt)-1l do begin
        ;; get the wav name
        restore,'data/specdata.sav'
        spawn,'cp data/mcmc/chains_'+datename+'/mcmc_chains_'+wavname[j]+'um.sav data/mcmc/mcmc_chains.sav'
        analyze_mcmc,discard=discardPoints
        spawn,'cp data/mcmc/param_unc/param_unc.txt data/mcmc/param_unc/param_unc_'+wavname[j]+'um.txt'
        firstwav = wavname[0]
        if keyword_set(docov) then begin
           analyze_cov,/psplot,discard=discardPoints
           spawn,'cp plots/mcmc/covar_plot.eps '+$
                 'plots/mcmc/individual_wavs/cov_plots_eps/'+datename+'/cov_plot_'+wavname[j]+'.eps'
        endif
     endfor
     gather_mcmc_radii
     plot_tim_ser,timebin=Npoints,/offtranserr,/psplot,/singleplot,deletePS=0,/showmcmc
     spawn,'cp plots/spec_t_series/tser_'+firstwav+'.eps plots/spec_t_series/tser_'+datename+'.eps'
;     try_mcmc
;     gather_mcmc_radii
;     plot_tim_ser,/showmcmc,timebin=Npoints,/offtranserr,/singleplot,/psplot,/pngcopy
     
;           spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
;                 'radius_vs_wavelength/straightdec23_leg00fit_freelimblin_mcmc_hypers_'+$
;                 'free_free_offtrans_err_009pt_modexp_kern_nyquist_2.txt'
;           spawn,'cp plots/mcmc/individual_wavs/histograms_png/* plots/mcmc/individual_wavs/histograms_png_dec23'
;           spawn,'cp plots/mcmc/individual_wavs/chain_plots_png/* plots/mcmc/individual_wavs/chain_plots_png_dec23'
;           spawn,'cp plots/mcmc/individual_wavs/cov_plots_png/* plots/mcmc/individual_wavs/cov_plots_png_dec23'
;           spawn,'cp -r radius_vs_wavelength/fit_data_mcmc/* radius_vs_wavelength/fit_data_mcmc_dec23/'
;           spawn,'cp plots/spec_t_series/tser_0.91.png plots/spec_t_series/all_t_series_dec23.png'
;           spawn,'cp plots/spec_t_series/tser_0.91.pdf plots/spec_t_series/all_t_series_dec23.pdf'
;           spawn,'cp data/mcmc/*um.sav data/mcmc/chains_dec23/'
     spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
           'radius_vs_wavelength/'+radcode+'_'+datename+'_leg00fit_freelimblin_mcmc_hypers_'+$
           'free_fixed_offtrans_err_009pt_modexp_kern.txt'
  endfor


end
