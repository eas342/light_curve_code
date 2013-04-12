pro analyze_mcmc,psplot=psplot,$
                 paramUpper=paramUpper,paramLower=paramLower,$
                 medparams=medparams
  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/mcmc/basic_mcmc'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
  endif


  ;; get the mcmc data
  restore,'data/mcmc/mcmc_chains.sav'
  ;;chainparams,lmfit,lmunct


  parnames = ['Rp/R*','b/R*','u1','u2','A/R*','A_0','A_1','A_2','A_3']
;  parfree =  [      1,     0,   1,   0,     0,    1,   1,    1,    0 ]
  freeInd = where(freep)
  nparams = n_elements(lmfit)
  nfree = total(freep)
  assert,n_elements(parnames),'=',nparams,'Warning # of parameter mismatch'

  ;; Returned parameters and uncertainties
  medparams = lmfit
  paramUpper = fltarr(nparams)
  paramLower = fltarr(nparams)

  !p.multi = [0,3,2]

  for i=0l,nfree-1l do begin
     pInd = freeInd[i];; parameter index
;     plothist,chainparams[pInd,*],xhist,yhist,bin=lmunct[pInd] * 8E,$
;              /nodata
;     yhist = histogram(chainparams[pInd,*],binsize=lmunct[pInd] *
;     0.1E,$
     if lmunct[pInd] EQ 0.0D then mybinsize = 0.001 else mybinsize = lmunct[pInd] * 0.5E
     yhist = histogram(chainparams[pInd,*],binsize=mybinsize,$
                       locations=xhistleft)
     xhist = xhistleft + lmunct[pInd] * 0.5E /2E
;     stop
     normalization = total(yhist)
     if n_elements(xhist) LE 1 then begin
        ;; IF the histogram is very sparse, just plot [0,0],[0,0]
        xhist = [-1,0,1]
        yhist = [0,1,0]
     endif
     plot,xhist,yhist,$
          xtitle=parnames[pInd],psym=10,$
          charsize=2,xmargin=[4.5,2.5],xticks=1,yticks=1,$
          yrange=[0,max(yhist)],ystyle=1,xrange=[min(xhist),max(xhist)],$
          xstyle=1,ymargin=[3.5,1.5],$
          xminor=5,xticklen=0.05,ytitle='Counts',ytickname=''
;          ymargin=[2,2]
     ;; Show the error bars from the Levenberg-Marquardt method
     oplot,lmfit[pInd] - [1E,1E] * lmunct[pInd],!y.crange,linestyle=2
     oplot,lmfit[pInd] + [1E,1E] * lmunct[pInd],!y.crange,linestyle=2
;     stop

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

  textcomment=string('Parameter','LM Fit','LM +/-','MCMC Fit','MCMC+','MCMC-',$
                    format='(A8,5(1x,A16))')
  forprint,parnames,lmfit,lmunct,medparams,paramUpper,paramLower,$
           format='(A8,5(1x,F16.5))',textout='data/mcmc/param_unc/param_unc.txt',$
           comment=textcomment

end
