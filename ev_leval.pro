function ev_leval,p,x=x,yin=y,yerr=yerr
;; Evaluates the likelihood function for a given set of
;; hyperparameters p and data (X,Y)
;; p[0] is the correlation function strength
;; p[1] is the correlation function time scale
;; p[2] is sigma - the error in the points


npts = n_elements(x)

;; Generate the correlation function
C = dblarr(npts,npts)
for i=0,npts -1l do begin
   for j=0l,npts-1l do begin
      Argument = -0.5D * ((x[i] - x[j])/p[1])^2
      if Argument LT -15D then C[i,j] = 0D else begin
         C[i,j] = p[0] * exp(-0.5D * ((x[i] - x[j])/p[1])^2)
      endelse
   endfor
endfor

;; add sigma to the correlation matrix
if n_elements(yerr) NE 0 then begin
   for i=0l,npts-1l do begin
      C[i,i] = C[i,i] + yerr[i]^2
   endfor
   
endif else begin
   for i=0l,npts-1l do begin
      C[i,i] = C[i,i] + p[2]^2
   endfor
endelse

;; Invert C
Cinv = invert(C)

;; Residuals
r = y

;; Factor out a sigma^2 and find the log of the determinant
if n_elements(yerr) NE 0 then begin
   sigFact = median(yerr)
endif else sigFact = p[2]

;; Find the log determinant from the method described by
;; http://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix/
CC = C
la_choldc,CC,/double
logdetermC = 2D * total(alog(diag_matrix(CC)))

;; 2 X Log Likelihood from Gibson et al. 2012, appendix A3
Likelihood = -(r ## Cinv ## transpose(r)) - logdetermC - double(npts) * 1.8378771D

;; minimize -L to maximize L
return,-Likelihood[0]

end
