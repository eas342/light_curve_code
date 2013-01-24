pro compile_spec,extraction2=extraction2
;; Compiles the spectra into a few simple arrays to look at the spectrophotometry

;Nwavbins = 35 ;; number of wavelength bins
Nwavbins = 9 ;; number of wavelength bins
SigRejCrit = 3 ;; number of sigma to reject when binning

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

;; Get the file names
if keyword_set(extraction2) then begin
   readcol,'file_lists/full_speclist.txt',filen,format='(A)'
endif else begin
   readcol,'file_lists/extraction1_full_list.txt',filen,format='(A)'
endelse
nfile = n_elements(filen)

;; Get the detector info
readcol,'data/detector_info.txt',skipline=1,$
       descrip,detectdata,format='(A,F)'
Gain = detectdata[0]
ReadN = detectdata[1]

;; Get the length of 
a = mrdfits(filen[0],0,header)
sizea = size(a)
Ngpts = sizea[1] ;; number of grid points

;; Make a grid for the wavelengths
DLam = fxpar(header,"CD1_1")
lamstart = fxpar(header,"CRVAL1")
lamgrid = (DLam * findgen(Ngpts) + lamstart)/1E4 ;; in microns

;; Make a grid for the flux, background and UT time
Nap = sizea[2] ;; number of apertures
flgrid = fltarr(Ngpts,Nap,nfile)
backgrid = fltarr(Ngpts,Nap,nfile)
utgrid = dblarr(nfile)

for i=0l,nfile-1l do begin
   ;; Read all files into the grid
   a2 = mrdfits(filen[i],0,header2,/silent)
   ;; Check that the number of points is corrrect
   sizea2 = size(a2)

   for j=0,Nap-1 do begin
      flgrid[*,j,i] = a2[*,j,0] * Gain
      backgrid[*,j,i] = a2[*,j,1] * Gain ;; multiply by gain
   endfor
   utgrid[i] = double(fxpar(header2,'MJD_OBS'))
endfor

;; Reset all zeros and negative flux values
badp = where(flgrid LE 0)
flgrid[badp] = !values.f_nan

;; Find the photon erros
ReadNarr = replicate(ReadN,Ngpts,Nap,Nfile)
ErrGrid = nansqrt( flgrid + backgrid + readnarr^2 )

;; Divide the two spectra
Divspec = flgrid[*,0,*] / flgrid[*,1,*]
fracE = nansqrt((ErrGrid[*,0,*]/flgrid[*,0,*])^2 + $
             (ErrGrid[*,1,*]/flgrid[*,1,*])^2 )
DivspecE = fracE * Divspec
SNR = Divspec / DivspecE

;; Do the wavelength binning for the divided spec
binfl = fltarr(Nwavbins,nfile)
binflE = fltarr(Nwavbins,nfile)
binGrid = (lamgrid[Ngpts-1] - lamgrid[0]) * findgen(Nwavbins)/float(Nwavbins-1) + $
          lamgrid[0]
binsizes = fltarr(Nwavbins) + (lamgrid[Ngpts-1] - lamgrid[0])/float(Nwavbins-1)

binind = fltarr(Nwavbins,Nap,nfile) ;; individual binned fluxes
binindE = fltarr(Nwavbins,Nap,nfile)

for i=0,nfile-1 do begin
   y = avg_series(lamgrid,Divspec[*,0,i],SNR[*,0,i],binGrid,binsizes,weighted=1,$
                  oreject=sigRejCrit,eArr=yerr,/silent,errIn=divSpecE[*,0,i])
   binfl[*,i] = y
   binflE[*,i] = yerr

   for k=0,Nap-1 do begin
      y2 = avg_series(lamgrid,flgrid[*,k,i],flgrid[*,k,i]/Errgrid[*,k,i],binGrid,$
                      binsizes,weighted=1,$
                      oreject=sigRejCrit,eArr=yerr2,/silent,errIn=Errgrid[*,k,i])
      binind[*,k,i] = y2
      binindE[*,k,i] = yerr2
   endfor
endfor


if keyword_set(psplot) then begin
   device,/close
   cgPS2PDF,plotnm
   ;; return the plotting to normal
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

; Save all data
save,flgrid,lamgrid,utgrid,bingrid,binfl,binflE,$
     ErrGrid,SNR,Divspec,DivspecE,backgrid,SNR,$
     Nwavbins,binsizes,binind,binindE,$
     filename='data/specdata.sav'
end
