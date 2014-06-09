FUNCTION sec_eclipse,ph,d,b,AoR,periodP,binP=binP,smooth=smooth

; ph is the phase
; d is the secondary eclipse depth
; b is the impact parameter
; u1 and u2 are limb darkening parameters
;aoR is the a over r star ratio, but I may add that later
; PeriodP is the OPTIONAL period of the planet for use with binP
; binP is the keyword that bins the light curve in order to convolve
; the exposure time window with the theoretical light curve
; smooth convolves the time series with a boxcar representing the
; exposure time window

negpt = where(d LT 0E)
if negpt NE [-1] then d[negpt] = -d[negpt]
f = quadlc(ph,sqrt(d),b,0.0E,0.0E,AoR,periodP,binP=binP,smooth=smooth)
assert,n_elements(d),'=',1,'Secondary eclipse only set up for a single radius at a time'
if negpt NE [-1] then f = 2E - f

return,f ; return the fluxes
end
