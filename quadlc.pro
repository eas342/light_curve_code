FUNCTION quadlc,ph,p,b,u1,u2,AoR,periodP,binP=binP,smooth=smooth
assert,p,'>=',0,"Warning! Rp/Rs is < 0"
; ph is the phase
; p is the planet size
; b is the impact parameter
; u1 and u2 are limb darkening parameters
;aoR is the a over r star ratio, but I may add that later
; PeriodP is the OPTIONAL period of the planet for use with binP
; binP is the keyword that bins the light curve in order to convolve
; the exposure time window with the theoretical light curve
; smooth convolves the time series with a boxcar representing the
; exposure time window

if keyword_set(binP) then begin
   expt = 29.426D/(60D*24D); the exposure time, days
   lightc=dblarr(n_elements(ph))
   for i=0,n_elements(ph)-1 do begin
      ; choose n points to average the flux
      npts=256
      exptph = expt/periodP; exposure time, converted to phase
      phbin = ph[i] + (dindgen(npts)/double(npts) - 0.5D)*exptph
      tbin=AoR*(sin(phbin*2D*!DPI))
      zbin=sqrt(tbin^2+b^2)           ; put in the impact parameter
      occultquad,zbin,u1,u2,p,flbin,flunibin
      lightc[i]=mean(flbin) ; return the mean flux for this data point
   endfor
endif else begin
   ;; convert between stellar radii and orbital phase for a large circular orbit
   t=AoR*(sin(ph*2D*!DPI))
   bp=b*cos(ph*2D*!DPI) ; project the impact parameter with the orbit
   z=sqrt(t^2+bp^2)              ; put in the projected impact parameter
   occultquad,z,u1,u2,p,lightc,luni ; calculate the light curve
   ;; luni is the uniform source light curve
endelse

if keyword_set(smooth) then lightc=boxcar(ph,lightc,periodP)


return,lightc ; return the fluxes
end
