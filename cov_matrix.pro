function cov_matrix,npts,x,theta0,theta1
;; Covariance matrix used by likelihood evaluation, simulated series
;; and plot_tim_ser

firstX = rebin(x,npts,npts)
secondX = transpose(firstX)

return,cov_kernel(firstX - secondX,theta0,theta1)
;Argument = -abs(firstX - secondX)^2 * theta1^2
;bigNeg = where(Argument LT -15D,complement=smallNeg)
;if bigNeg NE [-1] then C2[bigNeg] = 0D
;if smallNeg NE [-1] then C2[smallNeg] = theta0 * exp(Argument[smallNeg])/theta1^2

;seconds1 = systime(/seconds)
;
;C = fltarr(npts,npts)
;  for i=0,npts -1l do begin
;     for j=0l,npts-1l do begin
;        Argument = -abs((x[i] - x[j])) * theta1
;;        Argument = -((x[i] - x[j]) * theta1)^2
;        if Argument LT -15D then C[i,j] = 0D else begin
;;         C[i,j] = p[0]^2 * exp(Argument) * (1D + abs(x[i] - x[j])*p[1])
;           C[i,j] = theta0 * exp(Argument) / theta1
;;           C[i,j] = theta0 * exp(Argument)
;        endelse
;     endfor
;  endfor
;
;seconds2 = systime(/seconds)
;print,seconds2 - seconds1, ' seconds for first cov matrix'
;seconds3 = systime(/seconds)
;; Old code & timing

;seconds4 = systime(/seconds)
;print,seconds4 - seconds3,' seconds for second cov matrix'


end
