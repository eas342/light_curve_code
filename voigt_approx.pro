function voigt_approx,a,x,sigma

GHW = sigma * 1.17741 ;; Gaussian Half Width at Half Max
LHW = GHW * a ;; Lorentzian Half width at Half Max
sigV = 0.5346E * LHW + sqrt(0.2166E * LHW^2 + GHW^2)

d = (a - 1E)/(a + 1E)
cL = 0.68188 + 0.61293 * d - 0.18384 * d^2 - 0.11568 * d^3
cG = 0.32460 - 0.61825 * d + 0.17681 * d^2 + 0.12109 * d^3

;; Lorentzian Component
f = cL * 0.318310E * sigV/(X^2 + sigV^2)
ExpArgument =-0.693147E * X^2/sigV^2 
smallP = where(ExpArgument GT -18E)
if smallP NE [-1] then begin
   ;; Gaussian Component only applied to points that are within
   ;; e^(-18) the rest are assumed to be zero
   f[smallP] = f[smallP] + cG * 0.469719E * exp(expArgument[smallP])/sigV   
endif


return,f
end
