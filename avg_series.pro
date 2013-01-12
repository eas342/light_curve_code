function avg_series,dataUT,dataY,SNR,UTlist,explist,weighted=weighted,oreject=oreject,$
                    eArr=eArr,silent=silent,makestop=makestop,stdevArr=stdevArr,$
                    errIn=errIn
;; takes a series of UT times vs data values
;; and bins them into time bins
;; dataUT is the time of the data series
;; dataY is the value of the series
;; SNR is the Signal to ratio of the series (used for making weights)
;; UTlist are the start times of the bins (in hours)
;; explist are the exposure times of the bins (in hours)
;; the noweights keyword allows you to ignore the SNR and NOT do a
;; weighted average
;; oreject -- sets the number of sigma at which to reject an outlier. If
;; undefined, no outlier rejection occurs
;; eArr - an output vector that gives error estimates for the binned value
;; silent -- surpressses the output of the files
;; stdevArr -- an output vector that gives the standard deviation of
;;             the points within the bin

nser = n_elements(dataUT)
assert,nser,'=',n_elements(dataY),"Input array mismatch"
assert,nser,'=',n_elements(SNR),"Input array mismatch"

nUT = n_elements(UTlist)
assert,nUT,'=',n_elements(explist),"Bins array mismatch"
UTlistend = UTlist + explist    ;end times for exposures

binArr = dblarr(nUT)
eArr = dblarr(nUT)
stdevArr = dblarr(nUT)

if n_elements(silent) EQ 0 then vb = 1 else vb = 0 ;; vb for verbose

;; Calculate the weights in the manner of Taylor et al. 'Data Analysis'

if n_elements(weighted) EQ 0 then begin
   ws = SNR^2 ;; check if weighted is defined (default is true)
   if vb then print,'************* DOING WEIGHTED AVG **********'
endif else begin
   if weighted EQ 1 then begin
      ws = SNR^2 ;; check if weighted is true
      if vb then print,'************* DOING WEIGHTED AVG **********'
   endif else begin
      ws = dblarr(nser) + 1D 
      if vb then print,'************* DOING STRAIGHT UP AVG **********'
   endelse
endelse
if n_elements(oreject) NE 0 then begin
   if vb then print,'********** '+strtrim(oreject,1)+' Sigma Outliers will be Rejected *******'
endif

sigV = stddev(dataY,/nan)

for i=0l,nUT - 1l do begin
    ;;Sum the values based on the weights
    ind = where(dataUT GE UTlist[i] and dataUT LT UTlistend[i])

    ;; if outlier reject is on, then reject outliers
    if ind NE [-1] then begin
       if n_elements(oreject) EQ 0 then begin
          ;;       passp = lindgen(n_elements(ind)) ;; points that pass
          ;; make  sure that the SNR is defined
          passp = where(finite(SNR[ind]) EQ 1)
       endif else begin
          passp = where(abs(dataY[ind] - median(dataY[ind])) LT sigV*float(oreject) and finite(SNR[ind]) EQ 1)
       endelse
    endif else passp = [-1]

    if passp EQ [-1] then binArr[i] = !values.d_nan else begin

       binArr[i] = total(ws[ind[passp]] * dataY[ind[passp]],/nan)
       ;; the weights should only be added up if dataY is a real number
       ;; (not a NAN)
       binArr[i] = binArr[i] / total(ws[ind[passp]]*double(finite(dataY[ind[passp]])),/nan) 
       if n_elements(errIn) NE 0 then begin
          eArr[i] = sqrt(total( (ws[ind[passp]]*errIn[ind[passp]])^2,/nan ) )
          eArr[i] = eArr[i] / total(ws[ind[passp]]*double(finite(dataY[ind[passp]])),/nan)
       endif
       ;; error in the mean = stdev / sqrt(N)
       stdevArr[i] = stddev(dataY[ind[passp]],/nan) / sqrt(float(n_elements(passp)))
    endelse

endfor

return,  binArr

end
