function ev_mcmc,expr,X,Y,Yerr,start,chainL=chainL,parinfo=pi,maxp=maxp,$
                 hyperparams=hyperparams,noadjust=noadjust
;; Calculates the Light curve for a time series X and data Y using an
;; MCMC approach
;; Tried to use the same input as mpfit, X is independent variable,
;; expression is an expression of a function to be evaluated, start is
;; the starting variables pi is the parameter info (fixed vs. free
;; etc) and chainL is the chain length
;; maxp - maximum number of allowed points before exiting
;; hyperparams -- puts in hyper-parameters to describe errors and co-variance

  if n_elements(chainL) EQ 0 then chainL = 200l
  if n_elements(maxP) EQ 0 then maxP = 90000l
  nparams = n_elements(start)
  if n_elements(hyperparams) NE 0 then begin
     nhypers = n_elements(hyperparams.start)
     noHyperSwitch = 0
  endif else begin
     NoHyperSwitch = 1
     nhypers = 0l
  endelse 

  ;; update point where it shows a plot
  updatept = 100l

  fitparams = fltarr(nparams) ;; final best-fit parameters


  ;; array of parameter modifications uniform from -1 to +1
;  randUarray = 2E * randomu(0,nparams,maxP) - 1E
  randParray = randomn(0,nparams,maxP)
  if n_elements(hyperparams) NE 0 then begin
     randHarray = randomn(1,nhypers,maxP)
     randHarray[0,*] = randomn(4,maxP)
     ;; Here for the possibly exponentially distributed
     ;; parameter, let's try an exponential jump distribution
;     randHarray[1,*] = -1E * (alog(randomu(2,maxp)) + 1E)
     randHarray[1,*] = randomn(5,maxP)
     randHarray[2,*] = randomn(6,maxP)
  endif
  randKeeparr = randomu(100,maxP) ;; Keep threshholds

  chisQarray = fltarr(chainL) ;; chi-squared parameters

  freep = 1 - pi.fixed ;; free parameters
  nfree = total(freep)
  dof = float(n_elements(X) - nfree) ;; degrees of freedom

  ;; Start by doing a Levenberg-Marquardt minimum
  if total(pi[*].fixed) EQ nparams then begin
     result = start
     punct = fltarr(nparams)
  endif else result = mpfitexpr(expr,X,Y,Yerr,start,parinfo=pi,perr=punct)
  lmfit = result
  lmunct = punct

  chainparams = fltarr(nparams,chainL);; chain of parameters 
  chainparams[*,0] = result
;  chainparams[*,0] = [0.143,0.507,0.1,0,4.751,1.0,0,0,0]
  if n_elements(hyperparams) NE 0 then begin
     chainHypers = fltarr(nhypers,chainL)
     chainHypers[*,0] = hyperparams.start
  endif else chainHypers = [0E]
  
  modelY = expression_eval(expr,X,chainparams[*,0])
  if n_elements(hyperparams) NE 0 then begin
     chisQarray[0] = ev_leval(chainHypers[*,0],x=x,yin=(Y - ModelY),yerr=Yerr)
  endif else chisQarray[0] = total( ((Y - ModelY)/Yerr)^2 )

  j=1
  
  ;; use the rough covariance matrix to choose random bump sizes
;  jumps = rebin(punct,nparams,maxP) * randParray * 3E
  jumps = rebin(punct,nparams,maxP) * randParray * 0.8E
  if n_elements(hyperparams) NE 0 then begin
     hyperjumps = rebin(hyperparams.jumpsize,nhypers,maxP) * randHarray
     freeParams = where(pi.fixed NE 1)
     print,['Param 0'+strtrim(freeParams,1),'Hyper'+strtrim([0,1,2],1),'Chi-Sq']
     ;; Make the offset parameters jumps much larger when including hyper-parameters
     jumps[5:8,*] = jumps[5:8,*] * 40E
  endif
;  jumps = rebin(punct,nparams,maxP) * randParray * 1.5E

  for i=0l,maxP-1l do begin
;     modelY = expression_eval(expr,X,chainparams[*,j-1])
;     chisQ = total( ((Y - modelY)/Yerr)^2); / dof

     chainparams[*,j] = chainparams[*,j-1] + jumps[*,i]
     if n_elements(hyperparams) NE 0 then begin
        chainHypers[*,j] = chainHypers[*,j-1] + hyperjumps[*,i]
     endif

     newModel = expression_eval(expr,X,chainparams[*,j])

     if n_elements(hyperparams) NE 0 then begin
        ;; Hyper parameters less than zero should have -infinite
        ;; likelihood, so we will skip the ev_leval procedure for
        ;; those
        if where(chainHypers[*,j] LT hyperparams(*).limits[0]) EQ [-1] and $
           (chainHypers[1,j] LT hyperparams(1).limits[1] OR hyperparams(1).limited[1] EQ 0) then begin
           newchisQ = ev_leval(chainHypers[*,j],x=x,yin=(Y - newModel),yerr=Yerr)
        endif else newchisQ =exp(50D)
     endif else newchisQ = total( ((Y - newModel)/Yerr)^2 )

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

     if j mod 100 EQ 0 and j LE 902 and not keyword_set(noadjust) then begin
        ;; Adjust step sizes on hyper-parameters to be roughly 1 sigma
        for k=0,nhypers-1l do begin
           madhyper = stddev(chainHypers[k,0:j-1])
           if madhyper NE 0 then begin
              hyperjumps[k,*] = hyperjumps[k,*]/mad(hyperjumps[k,*]) * madhyper * 0.5
           endif
        endfor
        ;; Adjust the step size on all regular parameters to be
        ;; roughly 1 sigma
        for k=0l,nparams-1l do begin
           madreg = stddev(chainparams[k,0:j-1])
           if madreg NE 0 then begin
              jumps[k,*] = jumps[k,*]/mad(jumps[k,*]) * madreg * 0.5
           endif
        endfor
     endif

     if i mod updatept EQ updatept-1l GT 0 then begin
        ;; Show a plot every now and then to show progress
        ;; Shorten the chains for plotting
        fullchain = chainparams
        chainparams = chainparams[*,0l:(j-1l)]
        aRatio = float(j)/float(i) ;; acceptance ratio        
        if n_elements(hyperparams) NE 0 then begin
           fullhypers = chainhypers
           chainhypers = chainhypers[*,0l:(j-1l)]
        endif
        save,chainparams,lmfit,lmunct,freep,dof,chisQarray,aRatio,$
             chainhypers,$
             filename='data/mcmc/mcmc_chains.sav'

        if j GT 10 then chainplot,nohyper=noHyperSwitch,/showL

        wait,0.02
        ;; Restore the chain to its full length
        chainparams = fullchain
        if n_elements(hyperparams) NE 0 then chainhypers = fullhypers
     endif

  endfor

  ;; Show the Covariance values from mpfit
  oplot,result[5] - [1E,1E] * punct[5],!y.crange,linestyle=2
  oplot,result[5] + [1E,1E] * punct[5],!y.crange,linestyle=2
  aRatio = float(j)/float(i) ;; acceptance ratio

  ;; Shorten the chains in case the run was truncated prematurely
  chainparams = chainparams[*,0l:(j-1l)]
  if n_elements(hyperparams) NE 0 then chainhypers = chainhypers[*,0l:(j-1l)]
  save,chainparams,lmfit,lmunct,freep,dof,chisQarray,aRatio,$
       chainhypers,$
       filename='data/mcmc/mcmc_chains.sav'
  ;; save the chains, the mpfit best-fit values
;  stop

  return,fitparams
end
