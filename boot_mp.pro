function boot_mp,expr,x,y,yerr,start,parinfo=parinfo,perr=perr,quiet=quiet,mperr=mperr,$
            diagnose=diagnose,ntry=ntry
;; this wrapper to mpfitexpr uses the bootstrap method to calculate
;; errors in the parameters

if n_elements(ntry) EQ 0 then ntry=200

nparams = n_elements(start)
npoints = n_elements(x)
paramArr = dblarr(nparams,ntry)
perr = fltarr(nparams)

if n_elements(seed) EQ 0 then seed=0
chooseInd = floor(randomu(seed,npoints,ntry) * float(npoints))

for i=0l,ntry-1l do begin
   paramArr[*,i] = mpfitexpr(expr,x[chooseInd[*,i]],y[chooseInd[*,i]],yerr[chooseInd[*,i]],$
                                                     start,parinfo=parinfo,/quiet)
endfor
result = mpfitexpr(expr,x,y,yerr,start,parinfo=parinfo,perr=mperr,/quiet)

for i=0l,nparams-1l do begin
   perr[i] = stddev(paramArr[i,*])

   if not keyword_set(quiet) then begin
      histbinsize = mperr[i] * 0.4E
      yhist = histogram(paramArr[i,*],binsize=histbinsize,locations=xhistleft)
      xhist = xhistleft + histbinsize/2E
      dat = struct_arrays(create_struct('xhist',xhist,'yhist',yhist))
      edat = create_struct('vertlines',mperr[i] * [1E,-1E] + result[i])
      print,'Mpfit Err= ',mperr[i],' Bootstrap Err= ',perr[i]
      
      genplot,dat,edat
      if quit_caught() then return,0E
   endif
endfor


return,result

end
