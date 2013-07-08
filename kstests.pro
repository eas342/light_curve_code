pro kstests,totsets=totsets
;; Performs KS tests between the spectra for different nights to
;; confirm consistency

  if n_elements(totsets) EQ 0 then totsets=2

  ;; read in the first radius versus wavelength file
  radfile = choose_file(searchDir='radius_vs_wavelength',$
                           filetype='.txt')
  ;; Search the radius_vs_wavlength directory for .txt files

  readcol,radfile,wavl,wavlsize,rad,rade,skipline=1,format='(F,F,F)'

  nwavs = n_elements(wavl)
  avgwavl = wavl
  avgwavlsize = wavlsize
  totrad = fltarr(nwavs,totsets)
  totrade = fltarr(nwavs,totsets)
  avgrad = fltarr(nwavs)
  avgrade = fltarr(nwavs)

  totrad[*,0] = rad
  totrade[*,0] = rade
     
  for i=1l,totsets-1l do begin
     print,'Choose Additional file ',strtrim(i-1l,1)
     file2 = choose_file(searchDir='radius_vs_wavelength',filetype='.txt')
     readcol,file2,wavl2,wavl2size,rad2,rade2,skipline=1,format='(F,F,F)'
     if n_elements(wavl2) NE nwavs then begin
        print,'Different number of wavelength bins in files!'
        return
     endif else begin
        assert,total(wavl2-avgwavl),'=',0,'Different Wavelength!'
        assert,total(wavl2size-avgwavlsize),'=',0,'Different Wavelength Sizes!'
        totrad[*,i] = rad2
        totrade[*,i] = rade2
     endelse
  endfor

  ;; compare the first file to the rest
   kstwo,totrad[*,0],totrad[*,1:totsets-1],D,prob
   print,'KS D ',D
   print,'KS Probability ',prob

  ;; Do some weighted averages
   weights = 1E / totrade^2
   avgrad0 = total(weights[*,0] * totrad[*,0]) / total(weights[*,0])
   avgrade0 = 1E / sqrt(total(weights[*,0]))
   avgradR = total(weights[*,1:totsets-1] * totrad[*,1:totsets-1]) / total(weights[*,1:totsets-1])
   avgradeR = 1E / sqrt(total(weights[*,1:totsets-1]))

   print,'Mean First Distrib: ',avgrad0,' +/- ',avgrade0
   print,'Mean of remaining: ',avgradR, ' +/- ',avgRadeR

end
