pro state_parameters,reInitialize=reInitialize,psplot=psplot,$
                     timSerRange=timSerRange,seeingDiff=seeingDiff,$
                     secondary=secondary,differential=differential
;; Plots the time series and then also a lot of other parameters below
;; to see if flux changes can be attributed to the FWHM changes,
;; drifts, airmass etc.
;; reInitialize -- find both the spectral shifts and FWHMs
;; psplot -- save a postscript plot
;; timSerRange -- show a custom time series range
;; seeingDiff - show the difference in seeing between the two stars
;; secondary -- looks at secondary eclipse data
;; differential - find the differential parameters between the two stars

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/spec_t_series/state_param'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=14,decomposed=1,/color
  endif


  paramnames = ["Airmass","FWHM (px)","Relative Position (px)","Individual Flux",$
                "Spectral Shift (~px)","Voigt a"]
  nparams = n_elements(paramnames)

  !p.multi=[0,1,nparams+1]

  plot_tim_ser,/singlep,/skipReset,custyrange=TimSerRange,custsep=0.02,secondary=secondary

  ;; get the phase information
  restore,'data/timedata.sav'

  if keyword_set(reInitialize) then begin
     get_profile_widths
  endif else begin
     ;; get the spectral info
     get_profile_widths,/useSaved
  endelse

  ;; get the profile widths and locations
  restore,'data/prof_widths.sav'

  ;; get the spectral info
  restore,'data/specdata.sav'

  ;; get the spec list name for the observation
  restore,'data/used_date.sav'
  restore,'data/shift_data/shift_'+specfileListNamePrefix+'.sav'

  ;; Make a structure with all state parameters
  statePStruct = create_struct('JD',utgrid,'phase',tplot,'FluxRatio',binfl[0,*])

  for i=0l,nparams-1l do begin
     case paramnames[i] of 
        "Airmass": begin
           y = airmass
           ShowY2 = 0
           shorthand = 'airmass'
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
           shorthand = 'fwhm'
        end
        "Relative Position (px)": begin
           y = starLocations[*,0]
           y = y - median(y)
           y2 = starLocations[*,1]
           y2 = y2 - median(y2)
           showY2 = 1
           shorthand = 'position'
        end
        "Individual Flux": begin
           y = double(transpose(binind[0,0,*]))
           y = y/median(y)
           y2 = double(transpose(binind[0,1,*]))
           y2 = y2/median(y2)
           showY2 = 1
           shorthand = 'individualFlux'
        end
        "Spectral Shift (~px)": begin
           y = double(transpose(specShiftArr[0,*]))
           y2 = double(transpose(specShiftArr[1,*]))
           showY2 = 1
           shorthand = 'specshift'
        end
        "Voigt a": begin
           y = transpose(voigts[*,0])
           y2 = transpose(voigts[*,1])
           showY2 = 1
           shorthand = 'voigtdamp'
        end
        "Slit Model (fraction)": begin
           d1 = double(transpose(specShiftArr[0,*]))
           d2 = double(transpose(specShiftArr[1,*]))
           H = 20.49E /2E
           sigma1 = widths[*,0]
           sigma2 = widths[*,1]
           trans1 = gauss_slit(d1 - 0D,H,sigma1)
           trans2 = gauss_slit(d2 + 9D,H,sigma2)
           y = trans1 / trans2
           showY2 = 0
;           y = trans1
;           y2 = trans2
;           showY2 = 1
           shorthand = 'slitmodel'
        end
        else: y = tplot * 0E
     endcase
     
     if keyword_set(differential) and showY2 then begin
        showY2 = 0
        y = y2 - y
        myYtitle = 'Diff '+paramnames[i]
     endif else myYtitle = paramnames[i]
     
     if ShowY2 then begin 
        myYrange = threshold([y,y2],low=0.1,high=0.9,mult=0.4) 
        statePStruct = create_struct(statePStruct,shorthand+'1',y,shorthand+'2',y2)
     endif else begin
        myYrange = threshold(y,low=0.1,high=0.9,mult=0.4)
        statePStruct = create_struct(statePStruct,shorthand,y)
     endelse
     myXrange = !x.crange

     plot,tplot,y,ytitle=myYtitle,$
          xtitle="Orbital phase",$
          yrange=myYrange,$
          xrange=myXrange,xstyle=1,psym=4,ystyle=1
     if ShowY2 then begin
        ;; For some variables we show the background AND reference
        oplot,tplot,y2,color=mycol('blue'),psym=4
        al_legend,['Planet Host','Reference'],$
               /right,/bottom,linestyle=[0,0],$
               color=[!p.color,mycol('blue')],charsize=0.5
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

  write_csv,'data/state_parameters/full_parameters/'+specfileListNamePrefix+'.csv',$
            statePStruct,header=tag_names(statePstruct)
  save,statePStruct,$
       filename='data/state_parameters/full_parameters/'+specfileListNamePrefix+'.sav'
  
end
