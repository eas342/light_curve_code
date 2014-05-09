function gauss_slit,d,H,sigma
;; This script calculates the transmission through a Gaussian Slit
;; d is the distance from the center part of the slit
;; 2H is the slit width
;; sigma is the star's standard deviation

;; sqrt2 = 1.4142135624D
Flux = 0.5D * $
       (erf((d + H)/(1.4142135624D * sigma)) - erf((d - H)/(1.4142135624D * sigma)))

return,Flux

end
