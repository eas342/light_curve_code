pro state_parameters,skipInitialize=skipInitialize,psplot=psplot,$
                     timSerRange=timSerRange,seeingDiff=seeingDiff
;; Plots the time series and then also a lot of other parameters below
;; to see if flux changes can be attributed to the FWHM changes,
;; drifts, airmass etc.
;; skipInitialize -- skip the initialization process
;; psplot -- save a postscript plot
;; timSerRange -- show a custom time series range
;; seeingDiff - show the difference in seeing between the two stars

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/spec_t_series/state_param'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=14,decomposed=1,/color
  endif

  if not keyword_set(skipInitialize) then begin
     get_profile_widths
  endif

  paramnames = ["Airmass","FWHM (px)","Relative Position (px)","Raw Flux"]
  nparams = n_elements(paramnames)

  !p.multi=[0,1,nparams+1]

  plot_tim_ser,/singlep,/skipReset,custyrange=TimSerRange

  ;; get the phase information
  restore,'data/timedata.sav'

  ;; get the spectral info
  restore,'data/specdata.sav'

  ;; get the profile widths and locations
  restore,'data/prof_widths.sav'

  for i=0l,nparams-1l do begin
     case paramnames[i] of 
        "Airmass": begin
           y = airmass
           ShowY2 = 0
        end
        "FWHM (px)": begin
           if keyword_set(seeingDiff) then begin
              y1 = widths[*,0] * 2.35E
              y2 = widths[*,1] * 2.35E
              showY2 = 0
              y = y2 - y1
              badp = where(y LT -5)
              if badp NE [-1] then y[badp] = !values.f_nan
           endif else begin
              y = widths[*,0] * 2.35E
              y2 = widths[*,1] * 2.35E
              showY2 = 1
           endelse
        end
        "Relative Position (px)": begin
           y = starLocations[*,0]
           y = y - median(y)
           y2 = starLocations[*,1]
           y2 = y2 - median(y2)
           showY2 = 1
        end
        "Raw Flux": begin
           y = double(transpose(binind[0,0,*]))
           y2 = double(transpose(binind[0,1,*]))
           showY2 = 1
        end
        else: y = tplot * 0E
     endcase
     
     if ShowY2 then myYrange = threshold([y,y2],low=0.1,high=0.9) else begin
        myYrange = threshold(y,low=0.1,high=0.9)
     endelse
     myXrange = !x.crange

     plot,tplot,y,ytitle=paramnames[i],$
          xtitle="Orbital phase",$
          yrange=myYrange,$
          xrange=myXrange,xstyle=1,psym=4
     if ShowY2 then begin
        ;; For some variables we show the background AND reference
        oplot,tplot,y2,color=mycol('blue'),psym=4
        legend,['Planet Host','Reference'],$
               /right,/bottom,linestyle=[0,0],$
               color=[!p.color,mycol('blue')],charsize=0.3
     endif
  endfor

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 300% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

  !p.multi=0

end
