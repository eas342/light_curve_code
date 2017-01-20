pro manual_shift,npx,nindex=nindex,interval=interval,selfRef=selfRef,selfInd=selfInd
  ;; Manually shifts spectra and displays the two spectra and the
  ;; division
  ;; npx - the numbers of pixels to shift by
  ;; nindex - the filename index
  ;; interval - the interval, if set it displays a title with the
  ;;            interval, index and shift for the GUI
  ;; selfRef

  if n_elements(nindex) EQ 0 then nindex = 25;; default
  if n_elements(npx) EQ 0 then npx = 0; default

  ;; Make sure it compiles 
  restore,'data/specdata.sav'
  
  !p.multi=[0,1,2]

  if n_elements(interval) NE 0 then begin
     custTitle = 'Index: '+string(nindex,format='(I03)')+$
                 '  Shift: '+string(npx,format='(F6.3)')+$
                 '  Interval: '+string(interval,format='(F6.3)')

  endif

  if n_elements(selfind) EQ 0 then begin
     ;; If the use doesn't use a specific index, use the middle
     ;; one for a reference
     selfind = round(n_elements(utgrid)/2)+1l
  endif
  if keyword_set(selfref) then begin
     ;; If it's self referenced, choose a target 'master' spectrum
     FulltargSpec = flgrid[*,0,selfind] ;; the compile_spec procedure puts target first
     FullorigRefspec = flgrid[*,0,nindex]
     ratioRange=[0.6,1.6]
  endif else begin
     FulltargSpec = flgrid[*,0,nindex] ;; the compile_spec procedure puts target first
     FullorigRefspec = flgrid[*,1,nindex]
     ratioRange=[0,3]
  endelse

  targSpec = FulltargSpec / median(FulltargSpec)
  origRefspec = FullorigRefspec / median(FullorigRefspec)

  refspec = shift_interp(origRefspec,npx)

  plot,lamgrid,targSpec,yrange=[0,3]
  oplot,lamgrid,refspec,color=mycol('green')

  if keyword_set(selfref) then begin
     legNames = ['Target Ind '+strtrim(selfind,1),'Target Ind '+strtrim(nindex)]
  endif else begin
     legNames = ['Target','Reference']
  endelse
  al_legend,legNames,color=[!p.color,mycol('green')],linestyle=0
  
  divspec = targSpec/refspec
  plot,lamgrid,divspec,yrange=ratioRange,title=custTitle

  !p.multi=0

end
