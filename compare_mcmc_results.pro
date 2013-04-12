pro compare_mcmc_results,psplot=psplot
  
;  file1='data/mcmc/mcmc_chains_0.1times_mpfit_error_jumps_3000_pt_0.95um.sav'
  file1='data/mcmc/mcmc_chains_1.5times_mpfit_error_jumps_3000_pt_0.95um.sav'
  file2='data/mcmc/mcmc_chains_1.0times_mpfit_error_jumps_3000_pt_0.95um.sav'
  mcnames = ['1.5-Sigma Jumps','1-Sigma Jumps']
  
  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/mcmc/Rp_histo_compare'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
  endif

  ;; get the mcmc data
;  restore,'data/mcmc/mcmc_chains.sav'
  restore,file1
  ;;chainparams,lmfit,lmunct

  parnames = ['Rp/R*','b/R*','u1','u2','A/R*','A_0','A_1','A_2','A_3']
;  parfree =  [      1,     0,   1,   0,     0,    1,   1,    1,    0 ]
  freeInd = where(freep)
  nparams = n_elements(lmfit)
  nfree = total(freep)
  assert,n_elements(parnames),'=',nparams,'Warning # of parameter mismatch'

;  !p.multi = [0,3,2]

;  for i=0l,nfree-1l do begin
  for i=0l,0l do begin
     pInd = freeInd[i];; parameter index
;     plothist,chainparams[pInd,*],xhist,yhist,bin=lmunct[pInd] * 8E,$
;              /nodata
     yhist = histogram(chainparams[pInd,*],binsize=lmunct[pInd] * 0.1E,$
                       locations=xhist)
     yhist = yhist * 3.3
     normalization = total(yhist)
     plot,xhist,yhist,$
          xtitle=parnames[pInd],psym=10,$
          charsize=2,xmargin=[7,2.5],xticks=1,yticks=4,$
          yrange=[0,max(yhist)*1.8],ystyle=1,xrange=[min(xhist),max(xhist)],$
          xstyle=1,ymargin=[3.5,1],ytitle='Counts',$
          xminor=5,xticklen=0.05
;          ymargin=[2,2]
     ;; Show the error bars from the Levenberg-Marquardt method
     oplot,lmfit[pInd] - [1E,1E] * lmunct[pInd],!y.crange,linestyle=2
     oplot,lmfit[pInd] + [1E,1E] * lmunct[pInd],!y.crange,linestyle=2

     ;; get the mcmc data
;  restore,'data/mcmc/mcmc_chains_10times_mpfit_error_jumps_4000_pt_0.95um.sav'
     restore,file2
     yhist2 = histogram(chainparams[pInd,*],binsize=lmunct[pInd] * 0.1E,$
                        locations=xhist2)
     oplot,xhist2,yhist2,psym=10,color=mycol('red')
     legend,[mcnames[0],mcnames[1],'Mpfit Covariance'],linestyle=[0,0,2],$
            color=mycol(['black','red','black']),/right,/clear

;     stop
  endfor

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 450% '+plotprenm+'.pdf '+plotprenm+'.png'
     plotprenm = 'plots/mcmc/chain_examples'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=12, ysize=8,decomposed=1,/color
  endif

;  restore,'data/mcmc/mcmc_chains.sav'
  restore,file1

  plot,chainparams[0,*],yrange=[0.99 *min(xhist),1.01 *max(xhist)],$
       ytitle='Rp/R*',xtitle='Chain Index',ystyle=1
  restore,file2
;  restore,'data/mcmc/mcmc_chains_10times_mpfit_error_jumps_4000_pt_0.95um.sav'
  oplot,chainparams[0,*],color=mycol('red')
  
  oplot,!x.crange,lmfit[pInd] - [1E,1E] * lmunct[pInd],linestyle=2
  oplot,!x.crange,lmfit[pInd] + [1E,1E] * lmunct[pInd],linestyle=2

  legend,[mcnames[0],mcnames[1],'Mpfit Covariance'],linestyle=[0,0,2],$
         color=mycol(['black','red','black']),/right,/clear

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 450% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

end
