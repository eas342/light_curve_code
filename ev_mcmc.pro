function ev_mcmc,expr,X,Y,Yerr,start,chainL=chainL,parinfo=pi
;; Calculates the Light curve for a time series X and data Y using an
;; MCMC approach
;; Tried to use the same input as mpfit, X is independent variable,
;; expression is an expression of a function to be evaluated, start is
;; the starting variables, pi is the parameter info (fixed vs. free
;; etc) and chainL is the chain length


  if n_elements(chainL) EQ 0 then chainL = 50l
  if n_elements(maxP) EQ 0 then maxP = 10000l
  nparams = n_elements(start)


  fitparams = fltarr(nparams) ;; final best-fit parameters
  Keeps = bytarr(chainL) ;; array describing whether trial point was kept

  randUarray = randomu(0,maxP) ;; array of parameter modifications
  randKeeparr = randomu(1,maxP) ;; Keep threshholds

  nfree = nparams - total(pi.fixed)
  dof = float(n_elements(X) - nfree) ;; degrees of freedom

  ;; Start closer to minimum
  chainparams = fltarr(nparams,chainL);; chain of parameters  
  chainparams[*,0] = start
  chainparams[5,0] = 0.365
  j=1

  for i=0l,maxP-1l do begin
     modelY = expression_eval(expr,X,chainparams[*,j-1])
     chisQ = total( ((Y - modelY)/Yerr)^2) / dof

     chainparams[*,j] = chainparams[*,j-1]
     chainparams[5,j] = chainparams[5,j-1] + (randUarray[i] -0.5) * 0.05

     newModel = expression_eval(expr,X,chainparams[*,j])
     newchisQ = total( ((Y - newModel)/Yerr)^2 ) /dof

     DeltaChisQ = newChisQ - chisQ
     if DeltaChisQ GT 50E then Lratio = exp(-50E) else begin
        if DeltaChisQ LT -50E then Lratio = exp(50E) else begin
           Lratio = exp(-DeltaChisQ)
        endelse
     endelse

     if Lratio GT randKeeparr[i] then begin;; Probability of keeping the point
        j++
;     print,chisQ
;     print,'Offset   Delta Chi-Squared   Lratio    Keep?'
;     print,chainparams[5,i],DeltaChisQ,lratio,keeps[i]
        if j GE chainL then break
     endif else chainparams[*,j] = chainparams[*,j-1] ;;return to old point

     if i mod 300 EQ 299 and total(keeps) GT 0 then begin
        ;; Show a plot every now and then to show progress
        plot,chainparams(5,*),ystyle=16
        wait,0.02
     endif
  endfor

  plothist,chainparams(5,*),bin=0.001
  save,chainparams,filename='data/mcmc/offset_array2.sav'
  stop

  return,fitparams
end
