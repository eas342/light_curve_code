function ev_sigV,sigmaG,a
;; Converts the input Gaussian sigma and Voigt damping parameter into
;; the Voigt FWHM (using formulae from Liu et al. 2001)
;; constant is sqrt(2 * ln(2))
fG = 1.17741E * sigmaG
fL = a * fG
sigV = 0.5346E * fL + sqrt(0.2166E * fL^2 + fG^2)
return,sigV
end
