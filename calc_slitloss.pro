pro calc_slitloss,specshiftArr,usedate,widths,voigts,apkey,trans1,trans2
;; Calculates the slit loss for  the given parameters
;; returns the values in trans1 and trans2

  d1 = double(transpose(specShiftArr[0,*])) +1E
  d2 = double(transpose(specShiftArr[1,*])) +1E
  if strmatch(usedate,'*2012*') OR strmatch(usedate,'*2013*') then begin
     ;; Different plate scale for Aladdin vs Hawaii-2RG detectors
     H = 20.49/2E              
  endif else H = 30.5/2E
  sigma1 = widths[*,0]
  sigma2 = widths[*,1]
  vt1 = voigts[*,apkey[0]]
  vt2 = voigts[*,apkey[1]]
  trans1 = vslit_approx(d1,H,sigma1,median(vt1))
  trans2 = vslit_approx(d2,H,sigma2,median(vt2))



end
