pro compilespec
;; Compiles the spectra into a few simple arrays to look at the spectrophotometry


if keyword_set(psplot) then begin
   if keyword_set(divide) then begin
      plotnm = 'divided_stars.eps'
   endif else plotnm = 'two_corot_stars.eps'
   set_plot,'ps'
   !p.font=0
   device,encapsulated=1, /helvetica,$
          filename=plotnm
   device,xsize=14,ysize=10,decomposed=1,/color
endif

;; Overplot the two apertures in a spectrum
;a = mrdfits('../bigdog/bigdog0001.a.ms.d.fits',0,header)
;a = mrdfits('../bigdog/bigdog0002.a.ms.d.fits',0,header)
a = mrdfits('../bigdog/bigdog0200.a.ms.d.fits',0,header)
sizea = size(a)
Ngpts = sizea[1];; number of grid points

;; Get the wavelength grid info
DLam = fxpar(header,"CD1_1")
lamstart = fxpar(header,"CRVAL1")
lamgrid = (DLam * findgen(Ngpts) + lamstart)/1E4 ;; in microns
;lamgrid = (DLam * findgen(Ngpts) + lamstart) ;; in pixels

if keyword_set(divide) then begin
   sp1 = a[*,0,0]
   sp2 = a[*,1,0]
   zeropt = where(sp2 LE 0,complement=goodp)
   if zeropt NE [-1] then sp2[zeropt] = !values.f_nan
   y = sp1 / sp2
   ynam = "Flux Ratio"
   yran = [0,0.8]
endif else begin
   y = a[*,1,0]
   ynam = "Flux (DN)"
   yran = [0,0]
endelse
plot,lamgrid,y,$
;     xtitle="Wavelength (um)",$
     xtitle="Pixel",$
     ytitle=ynam,$
     yrange=yRan
;     xrange=[2.4,2.5]
if not keyword_set(divide) then begin
   oplot,lamgrid,a[*,0,0],color=mycol('red')
   legend,['Reference','Corot-1'],color=mycol(['black','red']),$
          psym=0,/right
endif

;plot,lamgrid,a[*,1,0]/a[*,0,0],$
;     xtitle="Wavelength (um)",$
;     ytitle="Flux Ratio"
;oplot,lamgrid,a[*,0,0],color=mycol('red')


if keyword_set(psplot) then begin
   device,/close
   cgPS2PDF,plotnm
   ;; return the plotting to normal
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

end
