pro ev_print_params,wavname,paramnames,result,punct,pi,k,$
                    skipsig=skipsig
;; Prints the labels, parameters and uncertainties
;; it is designed to scale the text widths to the appropriate
;; uncertainties
;; It only prints fitted parameters (not fixed ones)
  
;; Print the fit parameters options with string lengths
;; corresponding to the parameter uncertainties
;; skipsig -- skip printing the parameter uncertainties
nparams = n_elements(paramnames)
nwavs = n_elements(wavname)

stringLengths = intarr(nparams)
numberStarts = intarr(nparams)
print,'Wavl  ',format='(A,$)'
for j=0l,nparams-1l do begin
   if not pi[j].fixed then begin
      numberStarts[j] = round(alog10(min(punct[j,*])))
      stringLengths[j] = abs(numberStarts[j]) + 3
      print,paramnames[j],format='(A'+strtrim(stringLengths[j],1)+',$)'
      print,'  ',format='(A,$)'
   endif
endfor
print,''
for k=0l,nwavs-1l do begin
;; print the fit parameters and uncertainties
   for j=0l,nparams-1l do begin
      if not pi[j].fixed then begin
         if j EQ 0 then print,wavname[k],'  ',format='(A4,A2,$)'
         print,result[j,k],$
               format='(F'+strtrim(stringLengths[j],1)+$
               '.'+strtrim(max([0,-numberStarts[j]]),1)+$
               ',$)'
         print,'  ',format='(A,$)'
      endif
   endfor
   print,''
   if not keyword_set(skipsig) then begin
      print,'+/-      ',format='(A,$)'
      for j=0l,nparams-1l do begin
         if not pi[j].fixed then begin
            print,punct[j,k],$
                  format='(F'+strtrim(stringLengths[j],1)+$
                  '.'+strtrim(max([0,-numberStarts[j]]),1)+$
                  ',$)'
            print,'  ',format='(A,$)'
         endif
      endfor
      print,''
   endif
endfor

end
