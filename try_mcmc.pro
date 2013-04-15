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
  pi[5].fixed = 0 ;; free the offset
  pi[6].fixed = 0 ;; free the linear coefficient
;  pi[7].fixed = 0 ;; free the quadratic coefficient

  ;; Let the lim darkening and quadratic coefficient be negative
  pi[2].limited=[0,0]
  pi[7].limited=[0,0]

  ;; Model expression
  expr = 'quadlc(X,P[0],P[1],P[2],P[3],P[4])* (P[5] + X * P[6] + X^2 * P[7] + X^3 * P[8])'

  ;; Go through the cleaned time series
  cd,c=currentd
  fileopt = file_search(currentd+'/data/cleaned_tim_ser/*.txt')
  totfiles = n_elements(fileopt)
  for i=0l,n_elements(fileopt)-1l do begin
     trimst = strsplit(fileopt[i],'/',/extract)
     trimname = trimst(n_elements(trimst)-1l)
     namespl = strsplit(trimname,'_',/extract)
     wavname = namespl[n_elements(namespl)-2l]

     readcol,fileopt[i],$
             phase,fl,flerr,modelfl,resid,$
             format='(F,F,F,F,F)',skipline=1
     result = ev_mcmc(expr,phase,fl,flerr,start,parinfo=pi,chainL = 3000l,maxp=99000l)
;     result = ev_mcmc(expr,phase,fl,flerr,start,parinfo=pi,chainL = 3000l,maxp=90l)
     analyze_mcmc,/psplot
     ;; Save the chains
     spawn,'cp data/mcmc/mcmc_chains.sav data/mcmc/mcmc_chains_'+wavname+'.sav'
     ;; Save the histogram plot
     spawn,'cp plots/mcmc/basic_mcmc.pdf plots/mcmc/individual_wavs/histos_'+wavname+'.pdf'
     spawn,'cp plots/mcmc/basic_mcmc.png plots/mcmc/individual_wavs/histos_'+wavname+'.png'
     ;; Save the parameter uncertainties
     spawn,'cp data/mcmc/param_unc/param_unc.txt data/mcmc/param_unc/param_unc_'+wavname+'.txt'
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
