pro show_backg_profiles,choice,psplot=psplot,shortRow=shortRow
;; Shows how the backgroun spectrum varies as a function of spatial position
;; choice - which of the files in file_lists/images_for_profiles.txt
;;          to choose (a number 0,1 etc.)
;; psplot - set up a postscript plot that also makes PNG and PDF files
;; rowbyrow -- shows a narrower range of rows

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/background_profiles/back_profile'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=14, ysize=10,decomposed=1,/color
  endif

readcol,'file_lists/images_for_profiles.txt',imageList,format='(A)'

a = mrdfits(imageList[choice],0,origHeader)

NX = fxpar(origHeader,'NAXIS1')
NY = fxpar(origHeader,'NAXIS2')

nShow = 12l ;; number of rows to show
if keyword_set(shortRow) then begin
   RowStart = 55l
   linechoices = RowStart + lindgen(nShow) * 3l
endif else begin
   linechoices = round(findgen(nshow)/float(nshow) * NY)
endelse


showStart=195
showEnd = 220
separation = 0.1
if choice EQ 1 then begin
   myXrange=[195,220]
   myXstyle = 1
endif else begin
   myXrange=[195,220]
   myXstyle=1
endelse

nlabel = 5l ;; number of profiles to label
ncycle = round(float(nShow) / float(nlabel))


for i=0l,nShow-1l do begin
   xplot = lindgen(showEnd-showStart+1l)+showStart
   yplot = a[showStart:showEnd,linechoices[i]]
   ynorm = yplot/median(yplot)
   offset = separation * float(i)
   if i EQ 0l then begin
      plot,xplot,ynorm,$
           xtitle='Spectral Pixel',$
           ytitle='Normalized Background Flux - Offset',$
           yrange = [1E - separation * float(nShow+1),1E + separation],$
           /nodata,xrange=myXrange,xstyle=myXstyle,$
           xmargin=[10,10]
   endif
   oplot,xplot,ynorm - offset
   if i mod ncycle EQ 0 then begin
      ;; Label the row plotted
      ;; somewhere on the right of the graph
;      xstring = (!x.crange[1] - !x.crange[0])*0.9E + !x.crange[0]
      xstring=!x.crange[1]
      tabinv,xplot,xstring,Yind
      dataperYpix = (!y.crange[1] - !y.crange[0]) /((!y.window[1] - !y.window[0]) * !d.y_vsize)
      xyouts,xstring,ynorm[Yind] - offset - !D.Y_CH_SIZE * dataperYpix/2E,$
             'Row '+strtrim(linechoices[i]+1,1)
   endif

endfor

  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 250% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

end
