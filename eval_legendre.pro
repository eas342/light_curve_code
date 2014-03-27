function eval_legendre,X,C
;; This evaluates a Legendre polynomial centered on the middle
nterms = n_elements(C)
Y = dblarr(n_elements(X))
Xnorm = (2D * X - Max(X) - Min(X))/(Max(X) - Min(X) + 3D-16)
for i=0l,nterms-1l do begin
   Y = Y + Legendre(Xnorm,i) * C[i]
endfor

return,Y
end
