pro night_diff
;; Calculates the number of sigma away from the weighted average

readcol,'radius_vs_wavelength/avg_rp_rs.txt',wavl,wavlsize,rad,rade,$
        skipline=1,format='(F,F,F,F)'

radfile = choose_file(searchDir='radius_vs_wavelength',$
                      filetype='.txt')
readcol,radfile,wavl2,wavlsize2,rad2,rade2,skipline=1,format='(F,F,F)'

print,'wavl   #sigma from weighed av'
nwav = n_elements(wavl)
resultArray = fltarr(nwav,2)
resultArray[*,0] = wavl
resultArray[*,1] = (rad - rad2)/rade2
print,transpose(resultArray)

end
