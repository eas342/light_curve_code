pro find_shift_values,discreet=discreet,showfits=showfits,psplot=psplot,$
                      dec23=dec23
;; Straightens A Spectrum so that the X direction is wavelength
;; discreet - only move by discreet steps
;; showfits -- shows the fits to cross-correlations
;; dec23 -- looks at the Dec 23 data

if keyword_set(dec23) then begin
   img = mrdfits('../IRTF_UT2011Dec23/proc/bigdog/masterarc.fits',0,origHeade)
endif else img = mrdfits('../IRTF_UT2012Jan04/proc/bigdog/masterarc.fits',0,origHeader)

imgSize = size(img)
NX = imgSize[1]
NY = imgSize[2]
shiftArray = fltarr(NY)

recimg = fltarr(NX,NY) ;; rectified image

topspec = img[*,0]

;; set the up the PS plot
if keyword_set(psplot) then begin
   set_plot,'ps'
   !p.font=0
   if keyword_set(showfits) then plotprenm = 'plots/image_curvature/cross_cor' else begin
      plotprenm = 'plots/image_curvature/curv_func'
   endelse

   device,encapsulated=1, /helvetica,$
          filename=plotprenm+'.eps'
   device,xsize=14, ysize=10,decomposed=1,/color
endif

;; Remove NaNs (set to 0)
badp = where(finite(img) NE 1)
if badp NE [-1] then img[badp] = 0 ;; remove NaNs

;lagsize = 30l
lagsize = 30l
lagarray = lindgen(lagsize) - lagsize/2l
for i=0l,NY-1l do begin
   crosscor = c_correlate(topspec,img[*,i],lagarray)
   if keyword_set(discreet) then begin
      maxCross = max(crosscor,peakv)
   endif else begin
      fitExpr = 'P[0] * SinC(P[1] * (X -P[2])) + P[3] + P[4] *  EXP(-0.5E * ((X - P[2])/P[5])^2)'
;      fitExpr = 'P[0] * SinC(P[1] * (X -P[2])) * P[4] * EXP(-0.5E * ((X - P[2])/P[5])^2) + P[3]'
;      fitExpr = 'P[0] * SinC(P[1] * (X -P[2])) + P[3] + P[4] * X + P[5] * X^2
      startParams = [1,0.3,5,8,0,3]

;      fitExpr = 'P[0] * SinC(P[1] * (X -P[2]))'
;      startParams = [1,0.3,0]
      
      result=mpfitexpr(fitExpr,lagarray,crosscor,fltarr(lagsize)+0.1,$
                       startparams,/quiet)
;      if keyword_set(showfits) or i GE 251 then begin
;      if keyword_set(showfits) then begin
   endelse


   if keyword_set(discreet) then begin
      recimg[*,i] = shift(img[*,i],-lagarray[peakv]) 
      shiftarray[i] = -lagarray[peakv]
   endif else begin
      recimg[*,i] = shift_interp(img[*,i],-result[2])
      shiftarray[i] = -result[2]
   endelse

   if i GE 120 and keyword_set(showfits) then begin
      plot,lagarray,crosscor,ystyle=16,$
           xtitle='Lag (px)',ytitle='Cross Correlation'
      oplot,[-shiftarray[i],-shiftarray[i]],!y.crange,color=mycol('red')
      if not keyword_set(discreet) then begin
         ;; Show the fitted peak if using interpolation
         oplot,lagarray,expression_eval(fitExpr,lagarray,result),color=mycol('blue')
;         plot,img[*,i],/noerase,xstyle=4,ystyle=4
      endif
      stop
   endif

endfor



;; Show the shifts
Pos = findgen(NY) ;; position
if not keyword_set(showfits) then begin
   plot,Pos,shiftarray,xtitle='Row Number',ytitle='Shift Amount'
endif
;; Fit a polynomial to the shift curve

Npoly = 3
goodmask = where(Pos LT 140 OR $
                 Pos GT 200)
PolyTrend = poly_fit(pos[goodmask],shiftarray[goodmask],Npoly)
PolyMod = fltarr(NY)
for i=0l,Npoly do begin
   PolyMod = PolyMod + PolyTrend[i] * Pos^i
endfor
if not keyword_set(showfits) then begin
   oplot,Pos,PolyMod,color=mycol('green')
endif


;; Save the spectral shift parameters
forprint,Pos,PolyMod,comment='#Y Position (row)   Shift (px)',$
         textout='data/shift_data/shift_vals_from_arc.txt'


if keyword_set(psplot) then begin
   device, /close
   cgPS2PDF,plotprenm+'.eps'
   spawn,'convert -density 160% '+plotprenm+'.pdf '+plotprenm+'.png'
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

if keyword_set(dec23) then begin
   outFitsNm = '../IRTF_UT2011Dec23/proc/bigdog_rectified/masterarc.fits'
endif else outFitsNm = '../IRTF_UT2012Jan04/proc/bigdog_rectified/masterarc_cross_shift.fits'
writefits,outFitsNm,recimg,origHeader

end
