pro multiple_mcmc_nights

  ;; Make light curves, and do MCMC fitting for multiple nights
  ;; Make sure to empty the Cleaned_tim_ser folder first

  read,'Have you emptied the cleaned_time_ser folder?',junk

  read,'Have you emptied the param_unc folder?',junk

  for i=0l,2l do begin
     case i of
        0: compile_spec,/dec23,nwavbins=9,/nyquist
        1: compile_spec,nwavbins=9,/nyquist ;; Jan 04
        2: compile_spec,/dec29,nwavbins=9,/nyquist
     endcase
     if i LE 1 then Npoints = 100 else Npoints = 50
     plot_tim_ser,timebin=Npoints,/offtranserr
     try_mcmc
     gather_mcmc_radii
     plot_tim_ser,/showmcmc,timebin=Npoints,/offtranserr,/singleplot,/psplot,/pngcopy
     case i of
        0: begin
           spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
                 'radius_vs_wavelength/straightdec23_leg00fit_freelimblin_mcmc_hypers_'+$
                 'free_free_offtrans_err_009pt_modexp_kern_nyquist_2.txt'
           spawn,'cp plots/mcmc/individual_wavs/histograms_png/* plots/mcmc/individual_wavs/histograms_png_dec23'
           spawn,'cp plots/mcmc/individual_wavs/chain_plots_png/* plots/mcmc/individual_wavs/chain_plots_png_dec23'
           spawn,'cp plots/mcmc/individual_wavs/cov_plots_png/* plots/mcmc/individual_wavs/cov_plots_png_dec23'
           spawn,'cp -r radius_vs_wavelength/fit_data_mcmc/* radius_vs_wavelength/fit_data_mcmc_dec23/'
           spawn,'cp plots/spec_t_series/tser_0.91.png plots/spec_t_series/all_t_series_dec23.png'
           spawn,'cp plots/spec_t_series/tser_0.91.pdf plots/spec_t_series/all_t_series_dec23.pdf'
           spawn,'cp data/mcmc/*um.sav data/mcmc/chains_dec23/'
        end

        1: begin
           spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
                 'radius_vs_wavelength/straightjan04_leg00fit_freelimblin_mcmc_hypers_'+$
                 'free_free_off_trans_err_009pt_modexp_kern_nyquist_2.txt'
           spawn,'cp plots/mcmc/individual_wavs/histograms_png/* plots/mcmc/individual_wavs/histograms_png_jan04'
           spawn,'cp plots/mcmc/individual_wavs/chain_plots_png/* plots/mcmc/individual_wavs/chain_plots_png_jan04'
           spawn,'cp plots/mcmc/individual_wavs/cov_plots_png/* plots/mcmc/individual_wavs/cov_plots_png_jan04'
           spawn,'cp -r radius_vs_wavelength/fit_data_mcmc/* radius_vs_wavelength/fit_data_mcmc_jan04/'
           spawn,'cp plots/spec_t_series/tser_0.91.png plots/spec_t_series/all_t_series_jan04.png'
           spawn,'cp plots/spec_t_series/tser_0.91.pdf plots/spec_t_series/all_t_series_jan04.pdf'
           spawn,'cp data/mcmc/*um.sav data/mcmc/chains_jan04/'
        end
        2: begin
           spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
                 'radius_vs_wavelength/straightdec29_leg00fit_freelimblin_mcmc_hypers_'+$
                 'free_free_off_trans_err_009pt_modexp_kern_nyquist.txt'
           spawn,'cp plots/mcmc/individual_wavs/histograms_png/* plots/mcmc/individual_wavs/histograms_png_dec29'
           spawn,'cp plots/mcmc/individual_wavs/chain_plots_png/* plots/mcmc/individual_wavs/chain_plots_png_dec29'
           spawn,'cp plots/mcmc/individual_wavs/cov_plots_png/* plots/mcmc/individual_wavs/cov_plots_png_dec29'
           spawn,'cp -r radius_vs_wavelength/fit_data_mcmc/* radius_vs_wavelength/fit_data_mcmc_dec29/'
           spawn,'cp plots/spec_t_series/tser_0.91.png plots/spec_t_series/all_t_series_dec29.png'
           spawn,'cp plots/spec_t_series/tser_0.91.pdf plots/spec_t_series/all_t_series_dec29.pdf'
           spawn,'cp data/mcmc/*um.sav data/mcmc/chains_dec29/'
        end
     endcase
  endfor


end
