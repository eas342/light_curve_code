pro manual_shift,npx,nindex=nindex,interval=interval
  ;; Manually shifts spectra and displays the two spectra and the
  ;; division
  ;; npx - the numbers of pixels to shift by
  ;; nindex - the filename index
  ;; interval - the interval, if set it displays a title with the
  ;;            interval, index and shift for the GUI

  if n_elements(nindex) EQ 0 then nindex = 25;; default
  if n_elements(npx) EQ 0 then npx = 0; default

  restore,'data/specdata.sav'
  
  !p.multi=[0,1,2]

  if n_elements(interval) NE 0 then begin
     custTitle = 'Index: '+string(nindex,format='(I03)')+$
                 '  Shift: '+string(npx,format='(F6.3)')+$
                 '  Interval: '+string(interval,format='(F6.3)')

  endif

  FulltargSpec = flgrid[*,0,nindex] ;; the compile_spec procedure puts target first
  targSpec = FulltargSpec / median(FulltargSpec)
  FullorigRefspec = flgrid[*,1,nindex]
  origRefspec = FullorigRefspec / median(FullorigRefspec)

  refspec = shift_interp(origRefspec,npx)

  plot,lamgrid,targSpec,yrange=[0,3]
  oplot,lamgrid,refspec,color=mycol('green')

  al_legend,['Target','Reference'],color=[!p.color,mycol('green')],linestyle=0
  
  divspec = targSpec/refspec
  plot,lamgrid,divspec,yrange=[0,3],title=custTitle

  !p.multi=0

end
