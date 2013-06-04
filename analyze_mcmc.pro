pro analyze_mcmc,psplot=psplot,nohyper=nohyper,extend2lm=extend2lm
;; Takes the MCMC results and plots histograms for the parameters
;; psplot - generates postscript, png and pdf plots
;; nohyper - don't display hyper parameters - often used if
;;           hyper-parameters are fixed
;; extend2lm -- extend the X range to show all Levenberg-Marquardt results

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/mcmc/basic_mcmc'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif


  ;; get the mcmc data
  restore,'data/mcmc/mcmc_chains.sav'
;  restore,'data/mcmc/mcmc_chains_0.91um.sav'
;  restore,'data/mcmc/mcmc_chains_1.43um.sav'
  ;;chainparams,lmfit,lmunct

  parnames = ['Rp/R*','b/R*','u1','u2','A/R*','A_0','A_1','A_2','A_3']
;  parfree =  [      1,     0,   1,   0,     0,    1,   1,    1,    0 ]

  sizePchain = size(chainparams)
  nregular = sizePchain[1] ;; number of regular parameters
  ;; If there are hyper-parameters, plot those as well
  if n_elements(chainhypers) NE 0 AND not keyword_set(nohyper) then begin
     nparams = nregular + 2

     freep = [freep,1,1] ; make 2 hyperparametsr free
     parnames = [parnames,cgGreek('Theta')+['!D0!N','!D1!N']]
     fullchain = fltarr(nregular+2,sizePchain[2])
     fullchain[0:nregular-1,*] = chainparams
     fullchain[nregular:nparams-1l,*] = chainhypers[0:1,*]
     chainparams = fullchain
     medparams = [lmfit,0,0]
  endif else begin
     nparams = n_elements(lmfit)
     ;; Returned parameters and uncertainties
     medparams = lmfit
  endelse

  freeInd = where(freep)
  paramUpper = fltarr(nparams)
  paramLower = fltarr(nparams)

  nfree = total(freep)
  assert,n_elements(parnames),'=',nparams,'Warning # of parameter mismatch'


  !p.multi = [0,3,2]
  !X.Omargin = [4,2]

  for i=0l,nfree-1l do begin
     pInd = freeInd[i] ;; parameter index
     mybinsize = robust_sigma(chainparams[pInd,*]) * 0.5E

     if mybinsize EQ 0 then mybinsize=1 ;; can't use a zero bin size
     yhist = histogram(chainparams[pInd,*],binsize=mybinsize,$
                       locations=xhistleft)
     xhist = xhistleft + mybinsize /2E

     normalization = total(yhist)
     if n_elements(xhist) LE 1 then begin
        ;; IF the histogram is very sparse, just plot [0,0],[0,0]
        xhist = [-1,0,1]
        yhist = [0,1,0]
     endif
     if keyword_set(extend2lm) and pInd LE nregular-1l then begin
        ;; Show the LM results in case their many sigma away, give a
        ;; breathing room of (Max - Min)*0.05
        breathRoom = (max(xhist) - min(xhist))*0.05
        myXrange=[min([xhist,lmfit[pInd] - lmunct[pInd] - breathRoom]),$
                  max([xhist,lmfit[pInd] + lmunct[pInd] + breathRoom])]
     endif else myXrange=[min(xhist),max(xhist)]

     plot,xhist,yhist,$
          xtitle=parnames[pInd],psym=10,$
          charsize=2,xmargin=[5,5],$
          yrange=[0,max(yhist)],ystyle=1,xrange=myXrange,$
          xstyle=1,ymargin=[4,1.5],$
          xticklen=0.05,ytitle='Counts!C',$
          xtick_get=xtickvals,xtickformat='(A1)',$;; supress & save tick labels
          ytick_get=ytickvals,ytickformat='(A1)'
     twotick_labels,xtickvals,ytickvals

     ;; Show the error bars from the Levenberg-Marquardt method
     if pInd LE nregular-1l then begin
        oplot,lmfit[pInd] - [1E,1E] * lmunct[pInd],!y.crange,linestyle=2
        oplot,lmfit[pInd] + [1E,1E] * lmunct[pInd],!y.crange,linestyle=2
     endif

     ;; Show the 68% confidence limits from MCMC
     sortedp = sort(chainparams[pInd,*])
     chainL = n_elements(chainparams[pInd,*])
     thresh = 0.68
     upperp = sortedp[round((0.5E + thresh/2E) * chainL)]
     lowerp = sortedp[round((0.5E - thresh/2E) * chainL)]
     conflimits = chainparams[pInd,[lowerp,upperp]]
     oplot,replicate(conflimits[0],2),!y.crange,linestyle=3,color=mycol('blue')
     oplot,replicate(conflimits[1],2),!y.crange,linestyle=3,color=mycol('blue')
     
     medparams[pInd] = median(chainparams[pInd,*])
     paramLower[pInd] = medparams[pInd] - conflimits[0]
     paramUpper[pInd] = conflimits[1] - medparams[pInd]

  endfor

  print,'Acceptance Ratio = ',aRatio

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 450% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1

  endif
  !p.multi = 0
  !X.omargin = [0,0]

  textcomment=string('Parameter','LM Fit','LM +/-','MCMC Fit','MCMC+','MCMC-',$
                    format='(A8,5(1x,A16))')

  if n_elements(chainhypers) NE 0 AND not keyword_set(nohyper) then begin
     lmfull = [lmfit,0,0] ;; full Levenbergy-marquardt array (0 for the hyper-parameters)
     lmunctfull = [lmunct,0,0]
  endif else begin
     lmfull = lmfit
     lmunctfull = lmunct
  endelse

  forprint,parnames,lmfull,lmunctfull,medparams,$
           paramUpper,paramLower,$
           format='(A8,5(1x,F16.5))',textout='data/mcmc/param_unc/param_unc.txt',$
           comment=textcomment

  save,medparams,paramLower,paramUpper,filename='data/mcmc/param_unc/param_unc.sav'
  ;; Save the uncertainties and medians for use by analyze_cov

end
