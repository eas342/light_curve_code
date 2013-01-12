function nansqrt,num
  sqR = num * 0.0E
  fpt = where(finite(num) EQ 1 and num GT 0D,complement=nfpt)
  if fpt NE [-1] then begin
     sqR[fpt] = sqrt(num[fpt])
  endif
  if nfpt NE [-1] then begin
     sqR[nfpt]= dblarr(n_elements(nfpt))*!values.f_nan
  endif
  return,sqR
end


