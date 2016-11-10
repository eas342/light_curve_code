function planck_ratio,lamgrid,t1,t2,nonorm=nonorm
;; This function returns the ratio between two Planck functions
;; lamgrid - wavelength grid in angstroms
;; t1 - temperature of the first star
;;   
  angs = lamgrid * 1E4
  yout = planck(angs,t1)/planck(angs,t2)
  ;; renormalize
  if not keyword_set(nonorm) then begin
     yout = yout/median(yout) * 0.4E
  endif
  return,yout
end
