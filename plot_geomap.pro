pro plot_geomap,psplot=psplot
;; Plots the output of geomap which is compiled by the
;; compile_goemap_output script

  ;; Get the filename parameters
  readcol,'param_input/geomap_filenm.txt',filedescrip,filen,$
          format='(A,A)'

  ;; Read in the geomap data
  readcol,filen[0],/silent,$
           filenl,xshift,xrms,yshift,yrms,xrot,yrot,$
          format='(A,F,F,F,F,F,F)',skipline=1

  ;; get the UT times, JD times exptimes, etc
  readcol,filen[1],/silent,$
          filenl2,expt,time,date,JD,$
          format='(A,F,A,A,D)'

  ;; Get the transit times
  ;; get the transit times
  readcol,filen[2],skipline=1,$
          format='(A,A)',epoch,tepoch
  t1stcontact = date_conv(tepoch[0],'JULIAN')
  t4thcontact = date_conv(tepoch[1],'JULIAN')
  tmidtrans = (t1stcontact + t4thcontact)/2D

  trel = (JD - tmidtrans)*(24E*60E) ;; time relative to mid-transit in minutes

  if keyword_set(psplot) then begin
     plotnm = 'plots/centroid_plot.eps'
;     plotnm = 'plots/rotation_plot.eps'
     set_plot,'ps'
     !p.font=0
     device,encapsulated=1, /helvetica,$
            filename=plotnm
     device,xsize=14,ysize=10,decomposed=1,/color
  endif


;  cgplot,trel,yrot,$
;  cgplot,trel,yshift,$
  cgplot,trel,xshift,$
       xtitle='Time - Transit (min)',$
;       ytitle='Y Shift (px)'
       ytitle='X Shift (px)'
;       ytitle='Y Rot (deg)'

  ;; Print the RMS
  print,'Stdev Xshift (px) = ',stddev(xshift)
  print,'Stdev Yshift (px) = ',stddev(yshift)
  print,'Stdev Xrot (deg) = ',stddev(xrot)
  print,'Stdev Yrot (deg) = ',stddev(yrot)

  ;; Find the max/min X shifts and print the filenames
  maxSHVal = max(xshift,maxSHI)
  print,"Maximum Shift filename = ",filenl[maxSHI]
  minSHVal = min(xshift,minSHI)
  print,"Minimum Shift filename = ",filenl[minSHI]
;  plot,trel,xrot,$
;       xtitle='Time - Transit (min)',$
;       ytitle='X Rotation (deg)'

  if keyword_set(psplot) then begin
     device,/close
     cgPS2PDF,plotnm
     ;; return the plotting to normal
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif
  
end
