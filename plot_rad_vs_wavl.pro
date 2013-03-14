pro plot_rad_vs_wavl,psplot=psplot,showstarspec=showstarspec,$
                     nbins=nbins,custfile=custfile,$
                     showtheospec=showtheospec,choosefile=choosefile,$
                     totsets=totsets
;;psplot -- saves a postscript plot
;;showstarspec -- shows a star spectrum on the same plot
;;nbins -- number of points bo bin in Rp/R*
;;custfile -- chooses a custom radius vs wavelength file
;;showtheospec -- shows a theoretical transmission spectrum
;;choosefile -- choose a file from the radius_vs_wavelength directory
;;totsets -- optional keyword to specify the number of total sets of
;;           data to over-plot (eliminates the old additionalfile keyword)

  !x.margin = [13,14]
  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/rad_vs_wavl'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif

  ;; restore the star spectrum to show where telluric & star features
  ;; might be
  restore,'data/specdata.sav'

  ;; read in the radius versus wavelength file
  case 1 of
     keyword_set(custfile): radfile=custfile
     keyword_set(choosefile): begin
        radfile = choose_file(searchDir='radius_vs_wavelength',$
                              filetype='.txt')
        ;; Search the radius_vs_wavlength directory for .txt files
     end
     else: radfile='radius_vs_wavelength/radius_vs_wavl.txt'
  endcase

  readcol,radfile,wavl,rad,rade,skipline=1,format='(F,F,F)'

  if keyword_set(showstarspec) then ytempstyle=8 else ytempstyle=1

  if n_elements(wavl) EQ 1 then begin
     binsizes = binsizes
  endif else begin
     binsizes = fltarr(n_elements(wavl)) + wavl[1]-wavl[0]
  endelse


  if n_elements(nbins) NE 0 then begin
     ;; if asked to, bin the Rp/R*
     nnewwavl = ceil(float(n_elements(wavl))/float(nbins))
     newwavl = fltarr(nnewwavl)
     newrad = fltarr(nnewwavl)
     newrade = fltarr(nnewwavl)
     newbinsizes = fltarr(nnewwavl)
     for i=0l,nnewwavl-1l do begin
        subInd = lindgen(nbins)+i*nbins
        ;; check that you haven't gone over the number of points
        allowedsubpts = where(subInd LE n_elements(wavl) -1l)
        allowedpts = subInd[allowedsubpts]
        newwavl[i] = mean(wavl[allowedpts])
        weights = 1E/rade[allowedpts]^2
        newrad[i] = total(weights *rad[allowedpts])/total(weights)
        newbinsizes[i] = binsizes[i] * float(n_elements(allowedpts))
        newrade[i] = sqrt(1E /total(weights))
     endfor
     wavl = newwavl
     rad = newrad
     rade = newrade
     binsizes = newbinsizes
     ;; Save the binned wavelength file
     forprint,wavl,rad,rade,$
              textout='radius_vs_wavelength/binned_rp_rs.txt',$
              comment='# Wavelength(um)    Rp/R*    Rp/R* Error'
  endif else begin
  endelse

  plot,wavl,rad,$
       xtitle='Wavelength (um)',$
       ytitle='Rp/R*',$
;       ystyle=16,xstyle=1,$
       ystyle=ytempstyle,xstyle=1,xrange=[0.8,2.55],$
       yrange=[0.12,0.17],/nodata
;       yrange=[0.12,0.16],/nodata
;  oploterror,wavl,rad,rade
  wavlwidth = binsizes/2E
  oploterror,wavl,rad,wavlwidth,rade,psym=3
;                color=mycol('yellow') 
  

  ;; As in Gibson et al. 2012, show 3 scale heights around the
  ;; adopted Rp/R* from Jacob Bean et al. 2012
  scaleH = 0.00115E
;  scaleH = 0.00115E * 2E
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433],color=mycol(['red'])
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]+3E*scaleH,color=mycol(['red']),linestyle=2
  plots,[!x.crange[0],!x.crange[1]],[0.1433,0.1433]-3E*scaleH,color=mycol(['red']),linestyle=2
        
     
     
  prevXrange = !x.crange
  if keyword_set(showtheospec) then begin ;; show the theoretical transmission spectrum
     readcol,'../models/fortney_g10mps_2500K_isothermal.csv',theowav,theorad,$
             skipline=6,format='(F,F)'
     radToPlanet = 0.1037E ;; found by eye to get the approximate Corot bandpass correctly
     oplot,theowav,theorad *radToPlanet,color=mycol('blue')
  endif

  ;; if undefined, only show one set of data
  if n_elements(totsets) EQ 0 then totsets=1


  if totsets GT 1l then begin
     legnamearr = strarr(totsets)
     print,'Name for file '+radfile+' ?'
     tempnm=''
     read,tempnm,format='(A)'
     legnamearr[0l] = tempnm
     colorchoices = mycol(['black','orange','purple','blue'])
     ncolchoices = n_elements(colorchoices)
     colorarr = colorchoices[lindgen(totsets) mod ncolchoices]
  endif
     
  for i=2l,totsets do begin
     ;; If asked to, overplot another Rad/vs wavlength file
     print,'Choose Additional file ',strtrim(i-1l,1)
     file2 = choose_file(searchDir='radius_vs_wavelength',filetype='.txt')
     readcol,file2,wavl2,rad2,rade2,skipline=1,format='(F,F,F)'
     ;; find the bin width
     binsizes2 = fltarr(n_elements(wavl2)) + wavl2[1]-wavl2[0]
     wavlwidth2 = binsizes2/2E
     oploterror,wavl2,rad2,wavlwidth2,rade2,psym=3,$
                color=colorchoices[i-1l]
     print,'Legend Name for data from file '+file2+' ?'
     tempnm = ''
     read,tempnm,format='(A)'
     legnamearr[i-1l] = tempnm
  endfor
  if totsets GT 1l then begin
     legend,legnamearr,psym=1l+lonarr(totsets),color=colorarr
  endif

  if keyword_set(showstarspec) then begin
     ;; plot the source spectrum
     plot,lamgrid,flgrid(*,0,1),/noerase,xrange=prevXrange,ystyle=5,xstyle=1,$
          yrange=[-6E5,6E5],/nodata
     oplot,lamgrid,flgrid(*,0,1),color=mycol('blue')
     axis,yaxis=1,yrange=!y.crange,color=mycol('blue'),/ystyle,$
          ytitle='Raw Source Flux (DN)'
     
     !x.margin = [10.0,3.0]
  endif

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

end  
  
