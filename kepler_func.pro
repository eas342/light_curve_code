function kepler_func,x,p
;; A simple parametrized form of the Kepler light curve found by the
;; Eureqa program
  restore,'data/kepler_curves/phase_folded_kepler_all_kic1255.sav'
  f = interpol(kfluxS,phaseX,x)
  f = (f - 1D) * p + 1D
  return,f

end
