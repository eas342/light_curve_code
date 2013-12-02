pro rerun_multiple_mcmc,docov=docov,both=both,showLegResid=showLegResid
;; Takes the MCMC chain data and re-runs portions of it
  ;; Make light curves, and do MCMC fitting for multiple nights
  ;; Make sure to empty the Cleaned_tim_ser folder first

;; docov -- makes covariance plots for the parameters and saves the
;;          eps files
;; both -- do both photometry and spectroscopy
;; showLegResid -- shows the Legendre polynomial residuals (before
;;                 fancy GP + MCMC fitting)

  read,'Have you emptied the cleaned_time_ser folder?',junk

  read,'Have you emptied the param_unc folder?',junk

  discardpoints = 1000

  datename = ['dec23','jan04','dec29']

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
     cd,c=currentd
     fileopt = file_search(currentd+'/data/cleaned_tim_ser/*.txt')
     totfiles = n_elements(fileopt)
     for j=0l,n_elements(fileopt)-1l do begin
        ;; get the wav name
        restore,'data/specdata.sav'
        spawn,'cp data/mcmc/chains_'+datename[i]+'/mcmc_chains_'+wavname[j]+'um.sav data/mcmc/mcmc_chains.sav'
        analyze_mcmc,discard=discardPoints
        spawn,'cp data/mcmc/param_unc/param_unc.txt data/mcmc/param_unc/param_unc_'+wavname[j]+'um.txt'
        firstwav = wavname[0]
        if keyword_set(docov) then begin
           analyze_cov,/psplot,discard=discardPoints
           spawn,'cp plots/mcmc/covar_plot.eps '+$
                 'plots/mcmc/individual_wavs/cov_plots_eps/'+datename[i]+'/cov_plot_'+wavname[j]+'.eps'
        endif
     endfor
     gather_mcmc_radii
     plot_tim_ser,timebin=Npoints,/offtranserr,/psplot,/singleplot,deletePS=0,/showmcmc
     spawn,'cp plots/spec_t_series/tser_'+firstwav+'.eps plots/spec_t_series/tser_'+datename[i]+'.eps'

     spawn,'cp radius_vs_wavelength/mcmc_rad_vs_wavl.txt '+$
           'radius_vs_wavelength/'+radcode+'_'+datename[i]+'_leg00fit_freelimblin_mcmc_hypers_'+$
           'free_free_offtrans_err_009pt_matern_general.txt'
     spawn,'cp -r radius_vs_wavelength/fit_data_mcmc/* radius_vs_wavelength/fit_data_mcmc_'+datename[i]+'/'

     ;; Make a copy of the AutoCorrelation plots
     if keyword_set(showLegResid) then begin
        plot_tim_ser,timebin=Npoints,/offtranserr,/singleplot,legord=0,/fitcurve,/freelimblin
     endif
     analyze_resids,/psplot,/fast,/showkern
     for j=0l,n_elements(wavname)-1 do begin
        spawn,'cp plots/power_spectrum/acf_plot_'+wavname[j]+'.png plots/power_spectrum/'+$
              datename[i]+'_acf_png/'
        spawn,'cp plots/power_spectrum/acf_plot_'+wavname[j]+'.eps plots/power_spectrum/'+$
              datename[i]+'_acf_eps/'
     endfor

  endfor


end
