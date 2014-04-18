function fold_phase,t,secondary=secondary
;; Takes the phase and subtracts or adds the largest integer so
;; that it's sort of centered around a phase of 0
  if keyword_set(secondary) then foldtime = t mod 1.0D else begin
     foldtime = t - double(round(t))
  endelse
  return, foldtime
end
