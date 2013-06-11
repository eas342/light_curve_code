pro analyze_cov,psplot=psplot,nocontour=nocontour,nohyper=nohyper,$
                discard=discard
;; Make a plot of the co-variance of the MCMC parameters
;; psplot - make postscript & png files
;; nocontours -- skips the 68% and 95% contours
;; nohyper -- doesn't show the hyper-parameters
;; discard -- number of initial points to discard before using the chain

  ;; Use analyze_mcmc to get the parameter uncertainties
  analyze_mcmc,discard=discard
  restore,'data/mcmc/param_unc/param_unc.sav'
  unct = (paramLower + paramUpper)/2E
  fitpars = medparams

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
  ;;chainparams,fitpars,unct


  ;; medparams, paramLower and paramUpper
;  restore,'data/mcmc/mcmc_chains_1.43um.sav'


  parnames = ['R!Dp!N/R!D*','b/R!D*!N','u!D1!N','u!D2!N','a/R!D*!N','A!D0!N','A!D1!N','A!D2!N','A!D3!N']
  sizePchain = size(chainparams)
  nregular = sizePchain[1] ;; number of regular parameters (not including hypers)

  if n_elements(chainhypers) NE 0 AND not keyword_set(nohyper) then begin
     nparams = nregular + 2
     freep = [freep,1,1] ; make 2 hyperparametsr free
     parnames = [parnames,cgGreek('Theta')+['!D0!N','!D1!N']]
     fullchain = fltarr(nregular+2,sizePchain[2])
     fullchain[0:nregular-1,*] = chainparams
     fullchain[nregular:nparams-1l,*] = chainhypers[0:1,*]
     chainparams = fullchain
     medparams = [fitpars,0,0]
  endif else begin
     nparams = n_elements(fitpars)
     ;; Returned parameters and uncertainties
     medparams = fitpars
  endelse

  ;; Throw out the first discard points
  if n_elements(discard) NE 0 then begin
     truncchain = fltarr(nparams,sizePchain[2]-discard)
     truncchain = chainparams[*,discard:sizePchain[2]-1l]
     chainparams = truncchain
  endif
  

  freeInd = where(freep)

  nfree = total(freep)
  assert,n_elements(parnames),'=',nparams,'Warning # of parameter mismatch'

  !p.multi = [0,nfree-1,nfree-1]
  !Y.OMargin = [8,1.5]
  !X.OMargin = [16,1.5]

  for i=1l,nfree-1l do begin
     YpInd = freeInd[i];; parameter index for Y axis
     for j=0l,nfree-2l do begin
        XpInd = freeInd[j] ;; parameter index for X axis
;     if unct[pInd] EQ 0.0D then mybinsize = 0.001 else mybinsize = unct[pInd] * 0.5E
;     yhist = histogram(chainparams[pInd,*],binsize=mybinsize,$
;                       locations=xhistleft)
;     xhist = xhistleft + unct[pInd] * 0.5E /2E
;     stop
        ;; Set up the default plot parameters
        myXmargin = [0,0]
        myYmargin = [0,0]
        myXtickName = replicate(' ',3)
        myYtickName = replicate(' ',3)
        myXtitle=''
        myYtitle=''
        myXrange = [-3,3] * unct[xPind] + fitpars[xPind]
        myYrange = [-3,3] * unct[yPind] + fitpars[yPind]
        myNodata = 0
        myYstyle=1
        myXstyle=1

        if j GT i-1 then begin
           myNodata = 1 ;; Don't replicate the top right portion of the plot
           myYstyle=1+4
           myXstyle=1+4
        endif

        if j EQ nfree-1l then myXmargin[1] = 1.5 ;;Right side
        if i EQ 0 then myYmargin[1] = 1.5 ;; Top
;        if i EQ nfree-1l then begin  ;; Bottom
;           myXtickName = [string(myXrange[0],format='(F8.4)'),$
;                          string((myXrange[1]-myXrange[0])/2E + myXrange[0],format='(F8.4)'),$
;                                 ' ']
;           myXtickName = ['','',' ']
;           myXtitle = parnames[XpInd]
;        endif

        plot,chainparams[XpInd,*],chainparams[YpInd,*],psym=3,$
;             xtitle=myXtitle,ytitle=myYtitle,$
             charsize=2,xmargin=myXmargin,ymargin=myYmargin,$
             xtick_get=xtickvals,xtickformat='(A1)',$;; suppress & save tick labels
             ytick_get=ytickvals,ytickformat='(A1)',$
             ystyle=myYstyle,xstyle=myXstyle,$
             xticklen=0.05,$
             xrange=myXrange,yrange=myYrange,$
             nodata=myNodata
        ;; Show tick labels and axis for bottom row
        if i EQ nfree-1l then begin
           twotick_labels,xtickvals,ytickvals,/noY,xtitle=parnames[XpInd],/xmid,xorient=45
        endif
        if j EQ 0 then begin ;; Left side
           twotick_labels,xtickvals,ytickvals,/noX,ytitle=parnames[YpInd],/ymid
        endif
        
        if not keyword_set(nocontour) then begin
           
           ;; Now show a contour at 68% confidence
           binX = unct[XpInd]*0.3
           binY = unct[YpInd]*0.3
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

           if level2sig GE level1sig then begin
              mylevels=level1sig
           endif else mylevels = [level2sig,level1sig]
           contour,bin2,xcontour,ycontour,color=mycol('white'),/overplot,$
                   xrange=myXrange,yrange=myYrange,nodata=myNodata,$
                   levels=mylevels,thick=3
           contour,bin2,xcontour,ycontour,color=mycol('red'),/overplot,$
                   xrange=myXrange,yrange=myYrange,nodata=myNodata,$
                   levels=mylevels,thick=2
           
           ;; Show the 68% confidence intervals
;        sortedp = sort(chainparams[XpInd,*])
;        chainL = n_elements(chainparams[XpInd,*])
;        thresh = 0.68
;        upperp = sortedp[round((0.5E + thresh/2E) * chainL)]
;        lowerp = sortedp[round((0.5E - thresh/2E) * chainL)]
;        conflimits = chainparams[XpInd,[lowerp,upperp]]
;        oplot,replicate(conflimits[0],2),!y.crange,linestyle=3,color=mycol('blue')
;        oplot,replicate(conflimits[1],2),!y.crange,linestyle=3,color=mycol('blue')
        endif
        
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


end
