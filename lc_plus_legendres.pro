function lc_plus_legendres,X,P
;; This function evaluates the light curve plus a Legendre polynomial baseline
  ;; Ensure the Legendre argument is between -1 and 1
  normX = (2D * X - Max(X) - Min(X))/(Max(X) - Min(X) + 3D-16)
  ;; added in a small epsilon to keep Legendre function from
  ;; complaining (below -1 or above +1 due to numerical errors)
  return,quadlc(X,P[0],P[1],P[2],P[3],P[4]) * $
         (P[5] + $
          P[6] * Legendre(normX,1) + $
          P[7] * Legendre(normX,2) + $
          P[8] * Legendre(normX,3))
          
end
