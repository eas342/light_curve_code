function ev_mcmc,expr,X,Y,Yerr,start,chainL=chainL,parinfo=pi,maxp=maxp
;; Calculates the Light curve for a time series X and data Y using an
;; MCMC approach
;; Tried to use the same input as mpfit, X is independent variable,
;; expression is an expression of a function to be evaluated, start is
;; the starting variables, pi is the parameter info (fixed vs. free
;; etc) and chainL is the chain length


  if n_elements(chainL) EQ 0 then chainL = 200l
  if n_elements(maxP) EQ 0 then maxP = 90000l
  nparams = n_elements(start)

  ;; update point where it shows a plot
  updatept = 100l

  fitparams = fltarr(nparams) ;; final best-fit parameters

  ;; array of parameter modifications uniform from -1 to +1
;  randUarray = 2E * randomu(0,nparams,maxP) - 1E
  randParray = randomn(0,nparams,maxP)
  randKeeparr = randomu(100,maxP) ;; Keep threshholds

  chisQarray = fltarr(chainL) ;; chi-squared parameters

  freep = 1 - pi.fixed ;; free parameters
  nfree = total(freep)
  dof = float(n_elements(X) - nfree) ;; degrees of freedom

  ;; Start by doing a Levenberg-Marquardt minimum
  result = mpfitexpr(expr,X,Y,Yerr,start,parinfo=pi,perr=punct)

  chainparams = fltarr(nparams,chainL);; chain of parameters  
  chainparams[*,0] = result
  modelY = expression_eval(expr,X,chainparams[*,0])
  chisQarray[0] = total( ((Y - modelY)/Yerr)^2)
;  chainparams[*,0] = start

  j=1

  ;; use the rough covariance matrix to choose random bump sizes
;  jumps = rebin(punct,nparams,maxP) * randParray * 3E
  jumps = rebin(punct,nparams,maxP) * randParray * 0.1E
;  jumps = rebin(punct,nparams,maxP) * randParray * 1.5E

  for i=0l,maxP-1l do begin
;     modelY = expression_eval(expr,X,chainparams[*,j-1])
;     chisQ = total( ((Y - modelY)/Yerr)^2); / dof

     chainparams[*,j] = chainparams[*,j-1] + jumps[*,i]

     newModel = expression_eval(expr,X,chainparams[*,j])
     newchisQ = total( ((Y - newModel)/Yerr)^2 );; /dof

     DeltaChisQ = newChisQ - chisQarray[j-1]
     if DeltaChisQ GT 50E then Lratio = exp(-50E) else begin
        if DeltaChisQ LT -50E then Lratio = exp(50E) else begin
           Lratio = exp(-DeltaChisQ)
        endelse
     endelse
;     print,'Offset   Delta Chi-Squared   Lratio    Keep?'
;     print,chainparams[5,j],DeltaChisQ,lratio,(Lratio GT randKeeparr[i])

;     if i mod 10 eq 9 then stop

     if Lratio GT randKeeparr[i] then begin;; Probability of keeping the point
        chisQarray[j] = newchisQ
        j++
        if j GE chainL then break
     endif else chainparams[*,j] = chainparams[*,j-1] ;;return to old point

     if i mod updatept EQ updatept-1l GT 0 then begin
        ;; Show a plot every now and then to show progress
        plot,chainparams(5,0:j-1),ystyle=16
        wait,0.02

     endif

  endfor

  plothist,chainparams(5,*),bin=punct[5] * 3,$
           xtitle='Offset',ytitle='PDF'
  
  ;; Show the Covariance values from mpfit
  oplot,result[5] - [1E,1E] * punct[5],!y.crange,linestyle=2
  oplot,result[5] + [1E,1E] * punct[5],!y.crange,linestyle=2
  lmfit = result
  lmunct = punct
  aRatio = float(j)/float(i) ;; acceptance ratio

  ;; Shorten the chains in case the run was truncated prematurely
  chainparams = chainparams[*,0l:(j-1l)]
  save,chainparams,lmfit,lmunct,freep,dof,chisQarray,aRatio,$
       filename='data/mcmc/mcmc_chains.sav'
  ;; save the chains, the mpfit best-fit values
;  stop

  return,fitparams
end