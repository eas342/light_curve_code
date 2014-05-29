function vslit_approx,x,H,sigma,a

GHW = sigma * 1.17741 ;; Gaussian Half Width at Half Max
LHW = GHW * a ;; Lorentzian Half width at Half Max
sigV = 0.5346E * LHW + sqrt(0.2166E * LHW^2 + GHW^2)

d = (a - 1E)/(a + 1E)
cL = 0.68188 + 0.61293 * d - 0.18384 * d^2 - 0.11568 * d^3
cG = 0.32460 - 0.61825 * d + 0.17681 * d^2 + 0.12109 * d^3

sepN = (x - H)/sigV ;; negative separation point
sepP = (x + H)/sigV ;; positive separation point

;; Lorentzian Component
f = cL * 0.318310E * (atan(sepP) - atan(sepN))
;: Gaussian Component
f = f + cG * 0.5E * (erf(0.8325546D * sepP) - erf(0.8325546D * sepN))

return,f
end
