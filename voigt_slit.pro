function voigt_slit,d,H,w,ain
;; This script calculates the transmission of a Voigt profile with a Slit
;; d is the distance from the center part of the slit
;; 2H is the slit width
;; w is the stellar width
;; a is the Voigt Damping parameter

nds = n_elements(d)
nHs = n_elements(H)
nws = n_elements(w)
nas = n_elements(ain)
npts = max([nDs,nHs,nws,nas])

Flux = fltarr(npts)

;; make sure all are full arrays
if nds LT npts then d = replicate(d,npts)
if nHs LT npts then H = replicate(H,npts)
if nws LT npts then w = replicate(w,npts) * sqrt(2E) else w = w * sqrt(2E)
if nas LT npts then ain = replicate(ain,npts)

for i=0l,npts-1l do begin
   a = ain[i]
   save,a,filename='data/function_passing/voigtDamp.sav'
   Flux[i] = qsimp('Voigt_fixedA',(-H[i] -d[i])/w[i],(H[i] -d[i])/w[i],eps=1E-5)/1.7724539D
endfor
return,Flux

end
