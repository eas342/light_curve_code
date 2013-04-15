pro analyze_cov,psplot=psplot
  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/mcmc/covar_plot'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=17, ysize=12,decomposed=1,/color
  endif


  ;; get the mcmc data
  restore,'data/mcmc/mcmc_chains.sav'
;  restore,'data/mcmc/mcmc_chains_0.91um.sav'
  ;;chainparams,lmfit,lmunct

  parnames = ['R!Dp!N/R!D*','b/R!D*!N','u!D1!N','u!D2!N','a/R!D*!N','A!D0!N','A!D1!N','A!D2!N','A!D3!N']

  freeInd = where(freep)
  nparams = n_elements(lmfit)
  nfree = total(freep)
  assert,n_elements(parnames),'=',nparams,'Warning # of parameter mismatch'

  ;; Returned parameters and uncertainties
  medparams = lmfit
  paramUpper = fltarr(nparams)
  paramLower = fltarr(nparams)

  !p.multi = [0,nfree-1,nfree-1]
  !Y.OMargin = [4,1.5]
  !X.OMargin = [14,1.5]

  for i=1l,nfree-1l do begin
     YpInd = freeInd[i];; parameter index for Y axis
     for j=0l,nfree-2l do begin
        XpInd = freeInd[j] ;; parameter index for X axis
;     if lmunct[pInd] EQ 0.0D then mybinsize = 0.001 else mybinsize = lmunct[pInd] * 0.5E
;     yhist = histogram(chainparams[pInd,*],binsize=mybinsize,$
;                       locations=xhistleft)
;     xhist = xhistleft + lmunct[pInd] * 0.5E /2E
;     stop
        ;; Set up the default plot parameters
        myXmargin = [0,0]
        myYmargin = [0,0]
        myXtickName = replicate(' ',3)
        myYtickName = replicate(' ',3)
        myXtitle=''
        myYtitle=''
        myXrange = [-2.5,2.5] * lmunct[xPind] + lmfit[xPind]
        myYrange = [-2.5,2.5] * lmunct[yPind] + lmfit[yPind]
        myNodata = 0

        if j GT i-1 then myNodata = 1 ;; Don't replicate the top right portion of the plot

        if j EQ 0 then begin ;; Left side
;           myYtickName = [string(myYrange[0],format='(F8.4)'),$
;                          string((myYrange[1]-myYrange[0])/2E + myYrange[0],format='(F8.4)'),$
;                                 ' ']
           myYtickName = ['','',' ']
           myYtitle = parnames[YpInd]
        endif
        if j EQ nfree-1l then myXmargin[1] = 1.5 ;;Right side
        if i EQ 0 then myYmargin[1] = 1.5 ;; Top
        if i EQ nfree-1l then begin  ;; Bottom
;           myXtickName = [string(myXrange[0],format='(F8.4)'),$
;                          string((myXrange[1]-myXrange[0])/2E + myXrange[0],format='(F8.4)'),$
;                                 ' ']
           myXtickName = ['','',' ']
           myXtitle = parnames[XpInd]
        endif

        plot,chainparams[XpInd,*],chainparams[YpInd,*],psym=3,$
             xtitle=myXtitle,ytitle=myYtitle,$
             charsize=2,xmargin=myXmargin,ymargin=myYmargin,$
             xticks=2,yticks=2,xtickname=myXtickName,ytickname=myYtickName,$
             ystyle=1,xstyle=1,$
             xminor=4,xticklen=0.05,yminor=4,$
             xrange=myXrange,yrange=myYrange,$
             nodata=myNodata

        ;; Now show a contour at 68% confidence
        binX = lmunct[XpInd]*0.3
        binY = lmunct[YpInd]*0.3
        bin2=hist_2d(chainparams[XpInd,*],chainparams[YpInd,*],$
                    bin1=binX,bin2=binY,$
                    min1=myXrange[0],max1=myXrange[1],$
                    min2=myYrange[0],max2=myYrange[1])
        
        ;; Get dimensions of image
        bin2dims = size(bin2,/dim)
        Nxpoints = bin2dims[0]
        xcontour = findgen(Nxpoints) * binX + myXrange[0] + binX/2E ;; middle of bins
        NYpoints = bin2dims[1]
        ycontour = findgen(NYpoints) * binY + myYrange[0] + binY/2E ;; middle of bins

        ;; Find the 68% and 95% contours
        normbin = bin2 / float(total(bin2)) ;; PDF
        sortbin = sort(normbin)
        cumulativeT = total(normbin[sortbin],/cumulative)
        tabinv,cumulativeT,(1E -0.68),sig1pt
        tabinv,cumulativeT,(1E -0.95),sig2pt
        level1sig = bin2[sortbin[sig1pt]] ;; 68% confidence level
        level2sig = bin2[sortbin[sig2pt]] ;; 95% confidence level

        contour,bin2,xcontour,ycontour,color=mycol('white'),/overplot,$
                xrange=myXrange,yrange=myYrange,nodata=myNodata,$
                levels=[level2sig,level1sig],thick=3
        contour,bin2,xcontour,ycontour,color=mycol('red'),/overplot,$
                xrange=myXrange,yrange=myYrange,nodata=myNodata,$
                levels=[level2sig,level1sig],thick=2

        ;; Show the 68% confidence intervals
;        sortedp = sort(chainparams[XpInd,*])
;        chainL = n_elements(chainparams[XpInd,*])
;        thresh = 0.68
;        upperp = sortedp[round((0.5E + thresh/2E) * chainL)]
;        lowerp = sortedp[round((0.5E - thresh/2E) * chainL)]
;        conflimits = chainparams[XpInd,[lowerp,upperp]]
;        oplot,replicate(conflimits[0],2),!y.crange,linestyle=3,color=mycol('blue')
;        oplot,replicate(conflimits[1],2),!y.crange,linestyle=3,color=mycol('blue')


     endfor
  endfor

  !Y.OMargin = [0,0]
  !X.OMargin = [0,0]

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
