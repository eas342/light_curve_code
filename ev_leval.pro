function ev_leval,p,x=x,y=y
;; Evaluates the likelihood function for a given set of
;; hyperparameters p and data (X,Y)
;; p[0] is the correlation function strength
;; p[1] is the correlation function time scale
;; p[2] is sigma - the error in the points

npts = n_elements(x)

;; Generate the correlation function
C = fltarr(npts,npts)
;for i=0,npts -1l do begin
;   for j=0l,npts-1l do begin
;      C[i,j] = p[0] * exp(-0.5E * ((x[i] - x[j])/p[1])^2)
;   endfor
;endfor

;; Let's try fitting data to a line instead - here p[0] and
;;                                            p[1] will be the linear coefficients
ymodel = x * p[0] + p[1]

;; add sigma to the correlation matrix
for i=0l,npts-1l do begin
   C[i,i] = C[i,i] + p[2]
endfor

;; Invert C
Cinv = invert(C)

;; Residuals
r = (y - ymodel)/p[2]

;; 2 X Log Likelihood from Gibson et al. 2012, appendix A3
Likelihood = -(r ## Cinv ## transpose(r)) - determ(C) - float(npts) * 0.7982
;; minimize -L to maximize L
return,-Likelihood[0]

end
