pro ev_phot
;; Does some simple photometry on images

;; get the list of images
readcol,'file_lists/photometry_dec23_short.fits',filenm,format='(A)'
nfile = n_elements(filenm)

;; get the median_flat
flatimg = mrdfits('cals/median_flat.fits',0,hdr)
goodflpt = where(flatimg GT 0)
badflpt = where(flatimg LE 0)

firstHeadr = headfits(filenm[0])

Xlength = sxpar(firstHeadr,'NAXIS1')
Ylength = sxpar(firstHeadr,'NAXIS2')

;; Dark frame

for i=0l,0l do begin
;; get the image
   img = mrdfits(filenm[i],/dscale)

   ;; subtract the bias

   ;; zero out all negatives
   negpts = where(img LT 0)
   img[negpts] = 0

   img2 = fltarr(Xlength,Ylength)
   ;; flat field the image
   if goodflpt NE [-1] then begin
      img2[goodflpt] = float(img[goodflpt]) / flatimg[goodflpt] 
   endif
   if badflpt NE [-1] then begin
      img2[badflpt] = 0.0
   endif
   writefits,'cals/example_flattened'+string(i,format='(I04)')+'.fits',img2,hdr

endfor

end
