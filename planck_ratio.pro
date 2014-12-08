function planck_ratio,lamgrid,t1,t2
;; This function returns the ratio between two Planck functions
;; lamgrid - wavelength grid in angstroms
;; t1 - temperature of the first star
;;   
  angs = lamgrid * 1E4
  yout = planck(angs,t1)/planck(angs,t2)
  ;; renormalize
  yout = yout/median(yout) * 0.4E
  return,yout
end
