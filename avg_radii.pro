pro avg_radii,totsets=totsets
;;totsets -- optional keyword to specify the number of total sets of
;;           data to avg, default is 2

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

  ;; Do a weighted average
  for j=0l,nwavs-1l do begin
     weights = 1E / totrade[j,*]^2
     avgrad[j] = total(weights * totrad[j,*]) / total(weights)
     avgrade[j] = 1E / sqrt(total(weights))
  endfor


  forprint,avgwavl,avgwavlsize,avgrad,avgrade,$
           textout='radius_vs_wavelength/avg_rp_rs.txt',$
           comment='# Wavelength(um)  Bin size   Rp/R*    Rp/R* Error'

end  
  
