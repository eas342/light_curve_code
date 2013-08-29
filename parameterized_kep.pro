function parameterized_kep,x,p
;; A simple parametrized form of the Kepler light curve found by the
;; Eureqa program

  x2 = x - 0.05
  f = -3.2E * atan(0.00093E,x2^2 - x2 * atan(0.012E,1.1E + 26E * x2))
  f = f * 0.001E * p + 1E
  return,f

end
