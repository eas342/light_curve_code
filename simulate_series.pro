pro simulate_series,theta=theta,Npoints=Npoints,psplot=psplot,$
                    custYrange=custYrange
;; Makes a random time series with correlated noise
;; theta - is the set of hyperparameters that govern the correlated
;;         noise - theta[0] is the maximum correlation coefficient and
;;         theta[1] is the inverse timescale of the errors
;; Npoints - specifies the number of time points in the simulation,
;;           allowing the user to change the default
;; psplot -- generates a postscript plot
;; custYrange -- allows the user to specify a custom Y range

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/simulated_series/sim_series'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=10, ysize=7,decomposed=1,/color
  endif

  if n_elements(Npoints) EQ 0 then Npoints=600l

  X = findgen(Npoints)

  Nrealizations = 5l
  
  startseeds = lindgen(Nrealizations) + 1l
  colchoices = mycol(['red','purple','blue','orange','brown','black'])
  colorarray = colchoices[lindgen(Npoints) mod n_elements(colchoices)]
  stychoices = [0,2,3,5,7]
  stylearr = stychoices[lindgen(Npoints) mod n_elements(stychoices)]

  ;; Make a correlation matrix, C
  C = fltarr(Npoints,Npoints)

  ;; Try a squared exponential kernel from Danielski 2013
  ;; theta[0] is the maximum covariance
  ;; theta[1] is the inverse length scale
  if n_elements(theta) EQ 0 then theta = [0.2, 200]

  for i=0l,Npoints-1l do begin
     for j=0l,Npoints-1l do begin
        C[i,j] = theta[0] * exp(-0.5E *((X[i] - X[j])/theta[1])^2)
     endfor
  endfor

  ;; Add sigma to all diagonal elements
  sigma = 1E
  for i=0l,Npoints-1l do begin
     C[i,i] = C[i,i] + sigma
  endfor
  U = C ;; copy the correlation matrix
  ;; Find Cholesky decomposition - this will be multiplied by an
  ;;                               uncorrelated data array
  LA_CHOLDC,U

  ;; find the inverse of the correlation matrix
  if n_elements(custYrange) EQ 0 then custYrange = [-10,10]
  custYtitle=cgGreek('sigma')+'= '+string(sigma,format='(G5.2)')+',  '+cgGreek('theta')+'!D0!N= '+$
             string(Theta[0],format='(G5.2)')+',  '+cgGreek('theta')+'!D1!N= '+$
             string(Theta[1],format='(G7.2)')

  for j=0l,Nrealizations-1l do begin 

     RandomSer = randomn(startseeds[j],Npoints)
     ;; Multiply by matrix from correlation matrix
     Y = RandomSer ## U
     if j EQ 0l then begin
        plot,x,y,yrange=custYrange,$
             xtitle='Time (counts)',$
             ytitle='Flux',$
             title=custYtitle
     endif else begin
        oplot,x,y,linestyle=stylearr[j],color=colorarray[j]
     endelse

  endfor


  
;  steparray = lindgen(Npoints)
;  autoC = a_correlate(y,steparray)

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
