function eval_poly,X,C
;; This evaluates a polynomial
nterms = n_elements(C)
Y = dblarr(n_elements(X))
for i=0l,nterms-1l do begin
   Y = Y + X^i * C[i]
endfor

return,Y
end
