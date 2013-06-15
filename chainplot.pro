pro chainplot,psplot=psplot,nohyper=nohyper,extend2lm=extend2lm,$
              discard=discard
;; Takes the MCMC results and plots the parameters as a function of
;; time to check for convergence
;; psplot - generates postscript, png and pdf plots
;; nohyper - don't display hyper parameters - often used if
;;           hyper-parameters are fixed
;; extend2lm -- extend the X range to show all Levenberg-Marquardt results

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/mcmc/mcmc_chains'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif


  ;; get the mcmc data
  restore,'data/mcmc/mcmc_chains.sav'
;  restore,'data/mcmc/mcmc_chains_0.91um.sav'
;  restore,'data/mcmc/mcmc_chains_1.43um.sav'
  ;;chainparams,lmfit,lmunct

  parnames = ['R!Dp!N/R!D*!N','b/R!D*!N','u!D1!N','u!D2!N','A/R!D*!N','A!D0!N','A!D1!N','A!D2!N','A!D3!N']
;  parfree =  [      1,     0,   1,   0,     0,    1,   1,    1,    0 ]

  sizePchain = size(chainparams)
  nregular = sizePchain[1] ;; number of regular parameters
  ;; If there are hyper-parameters, plot those as well
  if n_elements(chainhypers) NE 0 AND not keyword_set(nohyper) then begin
     nparams = nregular + 2

     freep = [freep,1] ; make 2 hyperparametsr free
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


  !p.multi = [0,1,nfree]
  !X.Omargin = [11,1]
  !Y.omargin = [5,1]

  for i=0l,nfree-1l do begin
     pInd = freeInd[i] ;; parameter index

     if i LT nfree-1l then myXtitle='' else myXtitle='Steps'
     plot,chainparams[pInd,*],$
;          ytitle=parnames[pInd],psym=10,$
          charsize=2,xmargin=[5,5],$
          ystyle=16,$
;          yrange=[0,max(yhist)],ystyle=1,xrange=myXrange,$
          xstyle=1,ymargin=[0,0],$
          xticklen=0.05,xtitle=myXtitle,$
          xtick_get=xtickvals,xtickformat='(A1)',$;; supress & save tick labels
          ytick_get=ytickvals,ytickformat='(A1)'
     if i LT nfree-1l then noX=1 else noX=0
     twotick_labels,xtickvals,ytickvals,/ymid,ytitle=parnames[pInd],$
                    noX=noX

     ;; Show the error bars from the Levenberg-Marquardt method
     if pInd LE nregular-1l then begin
        oplot,!x.crange,lmfit[pInd] - [1E,1E] * lmunct[pInd],linestyle=2
        oplot,!x.crange,lmfit[pInd] + [1E,1E] * lmunct[pInd],linestyle=2
     endif

     ;; Show the 68% confidence limits from MCMC
     sortedp = sort(chainparams[pInd,*])
     if n_elements(sortedp) GT 15 then begin
        chainL = n_elements(chainparams[pInd,*])
        thresh = 0.68
        upperp = sortedp[round((0.5E + thresh/2E) * chainL)]
        lowerp = sortedp[round((0.5E - thresh/2E) * chainL)]
        conflimits = chainparams[pInd,[lowerp,upperp]]
        oplot,!x.crange,replicate(conflimits[0],2),linestyle=3,color=mycol('blue')
        oplot,!x.crange,replicate(conflimits[1],2),linestyle=3,color=mycol('blue')
     endif

     ;; Show where the points are discarded
     if n_elements(discard) NE 0 then begin
        oplot,replicate(discard,2),!y.crange,color=mycol('red')
     endif
     
  endfor


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
  !Y.omargin = [0,0]

end
