function fold_phase,t
;; Takes the phase and subtracts or adds the largest integer so
;; that it's sort of centered around a phase of 0
  foldtime = t - double(round(t))
  return, foldtime
end
