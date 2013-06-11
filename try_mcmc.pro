pro try_mcmc,psplot=psplot
;; Tests out my MCMC fit to a time series

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/mcmc/basic_mcmc'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
  endif


;  readcol,'data/cleaned_tim_ser/timeser_1.43um_.txt',$
;  readcol,'data/cleaned_tim_ser/timeser_0.91um_.txt',$

  ;; get the planet info
  readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
          skipline=1
  planetdat = create_struct('null','')
  for l=0l,n_elements(info)-1l do begin
     planetdat = create_struct(planetdat,info[l],data[l])
  endfor

  u1parm = 0.1E  ;; starting limb darkening parameters
  u2parm = 0.0E

  start=double([planetdat.p,planetdat.b_impact,u1parm,u2parm,$
                planetdat.a_o_rstar,1.0D,0D,0D,0D])
  pi = replicate({fixed:1, limited:[1,0], limits:[0.0E,0.0E]},9)
  ;; start with parameters fixed and then modify
  pi[0].fixed = 0 ;; free the planet radius
  pi[2].fixed = 0 ;; free the linear limb darkening
;  pi[3].fixed = 0 ;; free the quadratic limb darkening
  pi[5].fixed = 0 ;; free the offset
;  pi[6].fixed = 0 ;; free the linear coefficient
;  pi[7].fixed = 0 ;; free the second Legendre coefficient
;  pi[8].fixed = 0 ;; free the third Legendre coefficient

  ;; Let the limb darkening, quadratic and cubic coefficient be negative
  pi[2].limited=[0,0]
  pi[3].limited=[0,0]
  pi[6].limited=[0,0]
  pi[7].limited=[0,0]
  pi[8].limited=[0,0]

  ;; Model expression
;  expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* (P[5] + X * P[6] + X^2  * P[7] + X^3 * P[8])'
  expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* ( P[5] + '+$
         'Legendre((2E * X - Max(X) - Min(X))/(Max(X) - Min(X)),1) * P[6] + '+$
         'Legendre((2E * X - Max(X) - Min(X))/(Max(X) - Min(X)),2) * P[7] + '+$
         'Legendre((2E * X - Max(X) - Min(X))/(Max(X) - Min(X)),3) * P[8])'

  ;; set up the hyperparameters
  hyperpi = replicate({fixed:0, limited:[1,0], limits:[0.0E,0.0E],$
                      start:0E,jumpsize:0E},3)
;  hyperpi[*].start = [0.0028,1D-5,0.0024]
;  hyperpi[*].jumpsize = [0.0001,1E-6,0]
  hyperpi[*].start = [0.0005,0.2,0.002]
  hyperpi[*].jumpsize = [0.0001,0.05,0]

  ;; Go through the cleaned time series
  cd,c=currentd
  fileopt = file_search(currentd+'/data/cleaned_tim_ser/*.txt')
  totfiles = n_elements(fileopt)
  for i=0l,n_elements(fileopt)-1l do begin
;  for i=0l,0l do begin
     trimst = strsplit(fileopt[i],'/',/extract)
     trimname = trimst(n_elements(trimst)-1l)
     namespl = strsplit(trimname,'_',/extract)
     wavname = namespl[n_elements(namespl)-2l]

     readcol,fileopt[i],$
             phase,fl,flerr,modelfl,resid,$
             format='(F,F,F,F,F)',skipline=1

;     result = ev_mcmc(expr,phase,fl,flerr,start,parinfo=pi,chainL = 3000l,maxp=99000l)
     result = ev_mcmc(expr,phase,fl,flerr,start,parinfo=pi,chainL = 5000l,maxp=99000l,$
                      hyperparams=hyperpi)
;     result = ev_mcmc(expr,phase,fl,flerr,start,parinfo=pi,chainL = 200l,maxp=99000l,$
;                      hyperparams=hyperpi)
;     result = ev_mcmc(expr,phase,fl,flerr,start,parinfo=pi,chainL = 3000l,maxp=90l)
     analyze_mcmc,/psplot,discard=1000
     ;; Save the chains
     spawn,'cp data/mcmc/mcmc_chains.sav data/mcmc/mcmc_chains_'+wavname+'.sav'
     ;; Save the histogram plot
     spawn,'cp plots/mcmc/basic_mcmc.png plots/mcmc/individual_wavs/histograms_png/histos_'+wavname+'.png'
     spawn,'cp plots/mcmc/basic_mcmc.eps plots/mcmc/individual_wavs/histograms_eps/histos_'+wavname+'.eps'
     ;; Save the parameter uncertainties
     spawn,'cp data/mcmc/param_unc/param_unc.txt data/mcmc/param_unc/param_unc_'+wavname+'.txt'

     chainplot,/psplot,discard=1000
     ;; Save the chain plot
     spawn,'cp plots/mcmc/mcmc_chains.png plots/mcmc/individual_wavs/chain_plots_png/mcmc_chains_'+wavname+'.png'
     spawn,'cp plots/mcmc/mcmc_chains.eps plots/mcmc/individual_wavs/chain_plots_eps/mcmc_chains_'+wavname+'.eps'

     analyze_cov,/psplot,discard=1000
     ;; Save the parameter covariance plot
     spawn,'cp plots/mcmc/covar_plot.png plots/mcmc/individual_wavs/cov_plots_png/cov_plot_'+wavname+'.png'
     spawn,'cp plots/mcmc/covar_plot.eps plots/mcmc/individual_wavs/cov_plots_eps/cov_plot_'+wavname+'.eps'
  endfor



;  result = ev_mcmc(expr,phase,fl,flerr,start,parinfo=pi,chainl=600l)
;
  
  


  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif


end
