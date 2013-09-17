pro multiple_mcmc_nights,phot=phot,both=both
;; Make light curves, and do MCMC fitting for multiple nights
;; Make sure to empty the Cleaned_tim_ser folder first
;; phot - look at photometry data instead of spectroscopy
;; both - look at both photometry & spectroscopy

  read,'Have you emptied the cleaned_time_ser folder?',junk

  read,'Have you emptied the param_unc folder?',junk

  for i=0l,2l do begin
     case 1 of
        keyword_set(phot): begin
           case i of
              0: compile_phot,/dec23
              1: compile_phot
              2: compile_phot,/dec29
           endcase
           radCode = 'moris' ;; describing photometry data
        end
        keyword_set(both): begin
           case i of
              0: compile_both,/dec23
              1: compile_both
              2: compile_both,/dec29
           endcase
           radCode = 'both'
        end
        else: begin
           case i of
              0: compile_spec,/dec23,nwavbins=9
              1: compile_spec,nwavbins=9 ;; Jan 04
              2: compile_spec,/dec29,nwavbins=9
           endcase
           radCode = 'sspeX' ;; describing the data as spectra, (straightened)
        end
     endcase

     if i LE 1 then Npoints = 100 else Npoints = 50
     plot_tim_ser,timebin=Npoints,/offtranserr
     try_mcmc
     gather_mcmc_radii
     plot_tim_ser,/showmcmc,timebin=Npoints,/offtranserr,/singleplot,/psplot,/pngcopy
     case i of
        0: nightname = 'dec23'
        1: nightname = 'jan04'
        2: nightname = 'dec29'
     endcase
     
     ;; get the wavname so we know which time-ser file was created
     restore,'data/specdata.sav'

     spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
           'radius_vs_wavelength/'+radcode+'_'+nightname+'_leg00fit_freelimblin_mcmc_hypers_'+$
           'free_free_offtrans_err_009pt_modexp_kern.txt'
     spawn,'cp plots/mcmc/individual_wavs/histograms_png/* plots/mcmc/individual_wavs/histograms_png_'+nightname
     spawn,'cp plots/mcmc/individual_wavs/chain_plots_png/* plots/mcmc/individual_wavs/chain_plots_png_'+nightname
     spawn,'cp plots/mcmc/individual_wavs/cov_plots_png/* plots/mcmc/individual_wavs/cov_plots_png_'+nightname
     spawn,'cp -r radius_vs_wavelength/fit_data_mcmc/* radius_vs_wavelength/fit_data_mcmc_'+nightname+'/'
     spawn,'cp plots/spec_t_series/tser_'+wavname[0]+'.png plots/spec_t_series/all_t_series_'+nightname+'.png'
     spawn,'cp plots/spec_t_series/tser_'+wavname[0]+'.pdf plots/spec_t_series/all_t_series_'+nightname+'.pdf'
     spawn,'cp data/mcmc/*um.sav data/mcmc/chains_'+nightname+'/'
  endfor


end
