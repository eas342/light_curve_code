pro compile_multi_night,differential=differential,$
                        mnwavbins=mnwavbins,noremovelinear=noremovelinear,$
                        masktelluric=masktelluric,photonly=photonly
;; Puts together multiple nights of data into one time series
;; differential - do differential spectroscopy and re-read in the
;;                cleaned time series
;; mnwavbins - pass on to compile_spec the number of spectral
;;             wavelength bins
;; noremovelinear - by default, it removes the linear trend in each
;;                  data set before adding them together,
;;                  Noremovelinear will skip this step
;; masktelluric passes masktelluric onto compile_spec
;; photonly - only look at photometry

if n_elements(noremovelinear) EQ 0 then noremovelinear=0

;; get the list of speclists
  readcol,'file_lists/multi_night.txt',listFile,format='(A)',comment='#'

nNights = n_elements(listFile)

for i=0l,nNights-1l do begin
   choose_speclist,fchoice=listFile[i]

   ;; Get the spectral data, removing linear trends in the time series
;   compile_both,/readC,/removelinear,nwavbins=mnwavbins
   if n_elements(mnwavbins) EQ 0 then mnwavbins=9
   nwavbins= mnwavbins
   case 1 of 
      keyword_set(differential): begin
         compile_spec,/readC,removelinear=(1-noremovelinear),nwavbins=nwavbins,/specshift,$
                      masktelluric=masktelluric,/normalize,/quickread
      end
      keyword_set(photonly): begin
         compile_phot,/readC,removelinear=(1-noremovelinear)
      end
      else: begin
         compile_both,/readC,removelinear=(1-noremovelinear),nwavbins=nwavbins,/specshift,$
                      masktelluric=masktelluric,/normalize
      end
   endcase
   restore,'data/specdata.sav'

   if i EQ 0l then begin
      binflNew = binfl
      binflENew = binflE
      if n_elements(binind) NE 0 then binindNew = binind else binindNew = 0E
      if n_elements(binind) NE 0 then binindENew = binindE else binindENew = 0E
      utgridNew = utgrid
   endif else begin
      binflNew = transpose([transpose(binflNew),transpose(binfl)])
      binflENew = transpose([transpose(binflENew),transpose(binflE)])
      if n_elements(binind) NE 0 then binindNew = transpose([transpose(binindNew),transpose(binind)])
      if n_elements(binind) NE 0 then binindENew = transpose([transpose(binindENew),transpose(binindE)])
      utgridNew = [utgridNew,utgrid]
   endelse

endfor

binfl = binflNew ;; binned flux ratio
binflE = binflENew ;; standard error
binind = binindNew ;; individual binned fluxes for the two stars
binindE = binindENew
if mnwavbins EQ 1 then begin
   ;; If there is 1 wavelength bin, then the double transpose erases 1
   ;; dimension so it has to be put back
   nfileTot = n_elements(binind[0,*])
   binind = reform(binind,1,2,nfileTot)
   binindE = reform(binindE,1,2,nfileTot)
endif

utgrid = utGridNew

; Save all data
save,flgrid,lamgrid,bingrid,binfl,binflE,backdiv,$
     utgrid,itimegrid,wavname,$
     ErrGrid,SNR,Divspec,DivspecE,backgrid,$
     Nwavbins,binsizes,binind,binindE,filen,$
     airmass,altitude,backgrid,header,$
     filename='data/specdata.sav'
end
