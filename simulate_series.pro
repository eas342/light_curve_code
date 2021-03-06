pro simulate_series,theta=theta,Npoints=Npoints,psplot=psplot,$
                    custYrange=custYrange,sigma=sigma,$
                    Nrealizations=Nrealizations,autoc=autoc,$
                    modified=modified,residplot=residplot,$
                    showestimator=showestimator
;; Makes a random time series with correlated noise
;; theta - is the set of hyperparameters that govern the correlated
;;         noise - theta[0] is the maximum correlation coefficient and
;;         theta[1] is the inverse timescale of the errors
;; Npoints - specifies the number of time points in the simulation,
;;           allowing the user to change the default
;; psplot -- generates a postscript plot
;; custYrange -- allows the user to specify a custom Y range
;; Nrealizations -- specifies the number of realizations
;; autoc -- show the autocorrelation functions instead of the time series
;; modified -- multiplies the kernel by parameter theta1
;; residplot -- shows a plot much like the residual plots from data anlaysis
;; showestimator - show the autocovariance esimator, using
;;                 Wei's "Time Series Analysis"

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/simulated_series/sim_series'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=10, ysize=7,decomposed=1,/color
  endif

  if n_elements(Npoints) EQ 0 then Npoints=100l

  X = dindgen(Npoints)/double(Npoints-1l) * 0.2D - 0.1D

  if n_elements(Nrealizations) EQ 0 then Nrealizations= 5l
  
  startseeds = lindgen(Nrealizations) + 1l
  colchoices = mycol(['red','purple','blue','orange','brown','black'])
  colorarray = colchoices[lindgen(Nrealizations) mod n_elements(colchoices)]
  stychoices = [0,2,3,5,7]
  stylearr = stychoices[lindgen(Nrealizations) mod n_elements(stychoices)]

  ;; Try a squared exponential kernel from Danielski 2013
  ;; theta[0] is the maximum covariance
  ;; theta[1] is the inverse length scale
  if n_elements(theta) EQ 0 then theta = [0.0002, 5,0]

  ;; ensure that's it's double precision
  theta= double(theta)

  C = cov_matrix(Npoints,x,theta)


  ;; Add sigma to all diagonal elements
  if n_elements(sigma) EQ 0 then sigma = 0.002E
  sigma = double(sigma)
  for i=0l,Npoints-1l do begin
     C[i,i] = C[i,i] + sigma^2
  endfor

  ;; Try a Test Correlation matrix
;  C = [[1,0.6,0.3],[0.6,1,0.5],[0.3,0.5,1]]
;  Npoints = 3
  CC = C ;; copy the correlation matrix
  ;; Find Cholesky decomposition - this will be multiplied by an
  ;;                               uncorrelated data array

  LA_CHOLDC,CC,/double

  ;; Only keep the lower triangle & transpose
  U = fltarr(Npoints,Npoints)
  for i=0l,nPoints-1l do begin
     for j=0l,i do begin
        U[j,i] = CC[j,i]
     endfor
  endfor
  U = transpose(U)
  autoArray = fltarr(Npoints,Nrealizations)

  ntheta = n_elements(theta)
  thetaString = replicate(',  ',ntheta) + replicate(cgGreek('theta'),ntheta) + $
                '!D'+strtrim(indgen(ntheta),1)+'!N= '+string(Theta,format='(G7.2)')
  thetaString = string(thetaString,format='('+strtrim(ntheta,1)+'A)') ;; concatenate
  custYtitle=cgGreek('sigma')+'= '+string(sigma,format='(G6.2)')+thetaString

  ;; If asked to, show the auto-correlation functions instead of the
  ;; time series
  if keyword_set(autoC) then steparray = lindgen(Npoints)
  if keyword_set(psplot) then myCharsize=0.65 else myCharsize=1.0

  for j=0l,Nrealizations-1l do begin 

     RandomSer = randomn(startseeds[j],Npoints,/double)
     ;; Multiply by matrix from correlation matrix
     Y = RandomSer ## U


     case 1 of
        keyword_set(autoC): begin
           autoArray[*,j] = a_correlate(y,steparray,/cov)
           if j EQ 0l then begin
              if n_elements(custYrange) EQ 0 then begin
                 BothArrays = [transpose(C[0,*]),autoArray[*,j]]
                 custYrange = [min(BothArrays),max(bothArrays)]
              endif
              plot,steparray,autoArray[*,j],$
                   ytitle='Autocovariance',$
                   xtitle='Lag',yrange=custYrange,$
                   title=custYtitle,charsize=mycharsize,$
                   xmargin=[13,3]
           endif else begin
              oplot,steparray,autoArray[*,j],color=colorarray[j]
           endelse
           ;; Show the initial covariance
           if j EQ Nrealizations-1l then begin
              oplot,C[0,*],linestyle=0,colo=mycol('black'),thick=10
              oplot,C[0,*],linestyle=0,colo=mycol('blue'),thick=6

              ;; Find the average auto-correlation function
              avgAuto = fltarr(Npoints)
              for l=0l,Npoints-1l do begin
                 avgAuto[l] = mean(autoArray[l,*])
              endfor
              oplot,avgAuto,color=mycol('black'),linestyle=2,thick=6
              oplot,avgAuto,color=mycol('yellow'),linestyle=2,thick=3

              legendNames = ['Individual AC','Input Kernel','Ensemble Avg AC']
              legendColor = mycol(['purple' ,'blue'        ,'yellow'         ]);; foreground color
              legendBacks = mycol(['purple' ,'black'       ,'black'          ]);; background color
              legendStyle = [0              ,         0    ,           2     ] ;; linestyle
              legendBthic = [ 1             ,      10      ,           6     ] ; background thickness
              legendFthic = [ 1             ,      6       ,           3     ] ; foregroudn thickness

              if keyword_set(showestimator) then begin
                 AutoEstimator = auto_estimator(C)
                 oplot,AutoEstimator,color=mycol('black'),linestyle=2,thick=6
                 oplot,AutoEstimator,color=mycol('orange'),linestyle=2,thick=3
                 legendNames = [legendNames,'AC Estimator']
                 legendColor = [legendColor,mycol('orange')]
                 legendBacks = [legendBacks,mycol('black')]
                 legendBthic = [legenBthic,6]
                 legendFthic = [lgendFthic,3]
              endif
              

              al_legend,legendNames,color=legendBacks,$
                        thick=legendBthic,/right,linestyle=legendStyle,/clear,charsize=myCharsize
              al_legend,legendNames,color=legendColor,$
                        thick=legendFthic,/right,linestyle=legendStyle,charsize=myCharsize
           endif
        end
        keyword_set(residplot): begin
           if j EQ 0l then begin
              if n_elements(custYrange) EQ 0 then begin
                    custYrange = [-0.02,0.02]
              endif
              plot,x,y * 100E,xtitle='Orbital phase',$
                   ytitle='Relative Flux (%)',/nodata,$
                   yrange=custYrange * 100E,$
                   title=custYtitle,charsize=mycharsize
              dx = (x[1] - x[0])/2E
              oploterror,x,y*100E,$
                         fltarr(Npoints)+dx,(fltarr(Npoints)+sigma)*100E,$
                         psym=3,/nohat
           endif
        end
        else: begin
           if j EQ 0l then begin
              if n_elements(custYrange) EQ 0 then begin
                    custYrange = [-0.02,0.02]
              endif
              plot,x,y,yrange=custYrange,$
                   xtitle='Orbital Phase',$
                   ytitle='Flux',charsize=mycharsize,$
                   title=custYtitle;,xmargin=[7,7]
           endif else begin
              oplot,x,y,linestyle=stylearr[j],color=colorarray[j]
           endelse
        end
     endcase

  endfor

  ;; Save one of the series to a file
  forprint,x,y,$
     textout='data/simulated_series/simser.txt',$
     comment='# Time    Flux',/silent,$
           format='(D,D)'
  

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 300% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
;  plot,steparray[1l:Npoints-1l],autoC[1l:Npoints-1l],/xlog
;  stop

end
