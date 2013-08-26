pro compile_spec,extraction2=extraction2,sum=sum,nwavbins=nwavbins,$
                 dec23=dec23,dec29=dec29,nyquist=nyquist,extremeRange=extremeRange,$
                 maskwater=maskwater,custRange=custRange,widewatermask=widewatermask,$
                 cleanbyeye=cleanbyeye,specshift=specshift,starshift=starshift,$
                 custmask=custmask,molecbin=molecbin,trycurved=trycurved,$
                 matchgrid=matchgrid,readCurrent=readCurrent,skipBJD=skipBJD,$
                 masktelluric=masktelluric
;; Compiles the spectra into a few simple arrays to look at the spectrophotometry
;; extraction2 -- uses whatever spectra are in the data directory
;; sum -- uses the variance weighted (optimal) extraction by
;;            default, but with the sum keyword it will use the
;;            standard sum extraction from the IRAF reduction
;; nwavbins -- sets the number of wavelength bins to create
;; dec23 -- look at the dec23 data set (default is jan04)
;; readCurrent - reads from file_lists/current_speclist.txt
;; nyquist -- sample the wavelength bands at 2X bandwidth for Nyquist sampling
;; extremeRange -- chooses the minimum to maximum wavelength bins
;; custRange -- allows for a custom wavelength range
;; watermask - make a mask over the variable water
;; feature to take out some of the suspiciuos variability
;; widewatermask -- increases the size of the water mask
;; cleanbyeye -- a shortened file list where I've removed the
;;               bad spectra by eye
;; specshift -- use the shifting procedure where each specturm is
;;                shifted w/ cross-correlation
;; starshift -- allows integer shifts of the host star vs reference star
;; custmask -- creates a custom max over a specified wavelength rage
;;             such as [1.45,1.55]
;; molecbin -- bins the spectra for a given molecule instead of a grid
;; trycurved -- analyzes the spectra extracted from non-rectified
;;                dispersed images instead of straightened data
;; matchgrid -- matches the output grid to the original grid (no
;;              wavelength binning)
;; skipBJD -- skips the conversion from JD_UTC to BJD_TDB
;; masktelluric -- masks all telluric features

;Nwavbins = 35 ;; number of wavelength bins
;Nwavbins = 9 ;; number of wavelength bins
;Nwavbins = 20 ;; number of wavelength bins

if n_elements(Nwavbins) EQ 0 then Nwavbins = 9

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
if keyword_set(extraction1) then begin
   readcol,'file_lists/extraction1_full_list.txt',filen,format='(A)'
endif else begin
   ;; Otherwise use the current extraction
   case 1 of
      keyword_set(dec23):begin
         if keyword_set(cleanbyeye) then begin
            readcol,'file_lists/corot1_dec23cleaned_by_eye.txt',filen,format='(A)',$
                    stringskip='#'
         endif else readcol,'file_lists/dec23_straight_spec.txt',filen,format='(A)'
      end
      keyword_set(dec29): begin
         readcol,'file_lists/corot1_dec29fullspec.txt',filen,format='(A)'
      end
      else: readcol,'file_lists/corot1_jan_04_all_rectified.txt',filen,format='(A)'
   endcase
endelse

if keyword_set(trycurved) then begin
   if keyword_set(dec23) then begin
      readcol,'file_lists/corot1_dec23fullspec.txt',filen,format='(A)'
   endif else begin
      readcol,'file_lists/full_speclist.txt',filen,format='(A)'
   endelse
endif
if keyword_set(readCurrent) then begin
   readcol,'file_lists/current_speclist.txt',filen,format='(A)'
endif
nfile = n_elements(filen)

;; Get the detector info
readcol,'data/detector_info.txt',skipline=1,$
       descrip,detectdata,format='(A,F)'
Gain = detectdata[0]
ReadN = detectdata[1]

;;Get the coordinate info for the stars on the slit
readcol,'data/object_coordinates.txt',skipline=1,$
        objname,raHr,decDeg,format='(A,D,D)'
raDeg = raHr / 24D * 360D

;;Get the observatory
readcol,'data/observatory_info.txt',obscode,obsnm,$
        format='(A,A)',skipline=1

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
airmass = dblarr(nfile)
Altitude = dblarr(nfile,Nap) ;; Altitude of each star

if keyword_set(sum) then begin
   SpecKey = 0
endif else SpecKey = 1

for i=0l,nfile-1l do begin
   ;; Read all files into the grid
   a2 = mrdfits(filen[i],0,header2,/silent)
   ;; Check that the number of points is corrrect
   sizea2 = size(a2)

   for j=0,Nap-1 do begin
      flgrid[*,j,i] = a2[*,j,SpecKey] * Gain
      backgrid[*,j,i] = a2[*,j,2] * Gain ;; multiply by gain
   endfor
   utgrid[i] = double(fxpar(header2,'MJD_OBS'))
   utgrid[i] = utgrid[i] + double(fxpar(header2,'ITIME'))/(3600D * 24D)

   airmass[i] = double(fxpar(header2,'AIRMASS'))
   ;; find the differential airmass between the two stars
   utarray = dblarr(Nap) + utgrid[i]
   eq2hor,raDeg,decDeg,utarray,alt,az,obsname=obscode[0]
   for j=0,Nap-1 do begin
      altitude[i,j] = alt[j]
   endfor
endfor

;; Reset all zeros and negative flux values
badp = where(flgrid LE 0)
flgrid[badp] = !values.f_nan

;; Mask water if asked to
if keyword_set(maskwater) or keyword_set(widewatermask) or n_elements(custmask) NE 0 then begin
   if keyword_set(maskwater) then begin
      waterFirst=[1.34,1.37]
   endif else begin
      if keyword_set(widewatermask) then waterFirst=[1.33,1.42] else begin
         waterFirst = custmask
      endelse
;      waterFirst=[1.30,1.48]
;      waterFirst=[1.30,1.50]
   endelse
   watermask = where(lamgrid GT waterFirst[0] and lamgrid LE waterFirst[1])
   for i=0l,nfile-1l do begin
      for j=0,Nap-1 do begin
         flgrid[watermask,j,i] = !values.f_nan
      endfor
   endfor
endif

;; Mask water if asked to
if keyword_set(maskwater) then begin
   waterFirst=[1.34,1.37]
   watermask = where(lamgrid GT waterFirst[0] and lamgrid LE waterFirst[1])
   for i=0l,nfile-1l do begin
      for j=0,Nap-1 do begin
         flgrid[watermask,j,i] = !values.f_nan
      endfor
   endfor
endif

;; Mask all telluric features
if keyword_set(masktelluric) then begin
   Ntelluric = 4
   telluricW = [[1.1,1.15],[1.3,1.48],[1.78,1.95],[1.97,2.05]]
   for k=0,Ntelluric-1 do begin
      feature = telluricW[*,k]
      featuremask = where(lamgrid GT feature[0] and lamgrid LE feature[1])
      for i=0l,nfile-1l do begin
         for j=0,Nap-1 do begin
            flgrid[featuremask,j,i] = !values.f_nan
         endfor
      endfor
   endfor
endif

;; Find the photon erros
ReadNarr = replicate(ReadN,Ngpts,Nap,Nfile)
ErrGrid = nansqrt( flgrid + backgrid + readnarr^2 )

;; Star Shift, the default is to move the reference star by -1 pixel
if n_elements(starshift) EQ 0 then starshift = -1
xygrid = fltarr(Ngpts,nfile)
xygrid[*,*] = flgrid[*,1,*] ;; adjust the reference star
shiftedGrid = shift_interp(xygrid,starshift)
flgrid[*,1,*] = shiftedGrid

;; Shift arrays
if keyword_set(specshift) then begin
   ;; Align the stars within their bins
   
   if keyword_set(trystraight) then begin
      ;; since the straightend spec has some weird properties above
      ;; 2.47um because of shifting spectra, mask those out
      badp = where(lamgrid GT 2.47)
      flgrid[badp,*,*] = !values.f_nan
   endif

   for i=0l,Nap-1l do begin
      xyspec = fltarr(Ngpts,nfile)
      xyspec[*,*] = flgrid[*,i,*]
      ShiftedGrid1 = find_shifts(xyspec,/cutEnds)
      flgrid[*,i,*] = ShiftedGrid1
   endfor

   ;; Align the stars with each other
   medspec0 = fltarr(Ngpts)
   medspec1 = fltarr(Ngpts)
   for i=0l,Ngpts-1l do begin
      medspec0[i] = median(flgrid[*,0,*])
      medspec1[i] = median(flgrid[*,1,*])
   endfor
   lagarray = lindgen(20l) - 10l
   badp = where(finite(medspec0) EQ 0)
   if badp NE [-1] then medspec0[badp] = 0.0E
   badp = where(finite(medspec1) EQ 0)
   if badp NE [-1] then medspec1[badp] = 0.0E
   CrossC = c_correlate(medspec0,medspec1,lagarray)
   PolyFit = poly_fit(lagarray,crossC,2)
   ShiftStars = polyFit[1]/(-2E * polyFit[2])
   xyspec = fltarr(Ngpts,nfile)
   xyspec[*,*] = flgrid[*,0,*]
   shiftedGrid = shift_interp(xyspec,ShiftStars)
   flgrid[*,0,*] = shiftedGrid
endif

;; Divide the two spectra
Divspec = flgrid[*,0,*] / flgrid[*,1,*]
fracE = nansqrt((ErrGrid[*,0,*]/flgrid[*,0,*])^2 + $
             (ErrGrid[*,1,*]/flgrid[*,1,*])^2 )
DivspecE = fracE * Divspec
SNR = Divspec / DivspecE

;; Do the wavelength binning for the divided spec
if keyword_set(extremeRange) then begin
   StartWav = lamgrid[0]
   EndWav = lamgrid[Ngpts-1]
endif else begin
   if n_elements(custRange) EQ 0 then begin
      StartWav = 0.82E
      EndWav = 2.40E ;; wavelength range that has useful scientific information
   endif else begin
      StartWav = custRange[0]
      EndWav = custRange[1]
   endelse
endelse
;; These are the starts of the bins (not the middles)
if keyword_set(molecbin) then begin
   Nwavbins = 2 ;; in the molecule and outside the molecule
   molecule='C2H2'
   readcol,'molec_inputs/C2H2_bin_locations.txt',molecStarts,molecEnds,$
           skipline=1,format='(F,F)'
   nmolecWavs = n_elements(molecStarts)
   for i=0l,nmolecWavs-1l do begin ;; collect wavelengths
      molecP = where(lamgrid GT molecStarts[i] and lamgrid LE molecEnds[i])
      if molecP NE [-1] then begin
         if n_elements(molecInd) EQ 0 then molecInd = molecP else begin
            molecInd = [molecInd,molecP]
         endelse
      endif
   endfor
   ;; Find the complement
   byteTest = bytarr(Ngpts)
   byteTest[molecInd] = 1
   outmolec = where(byteTest EQ 0)

   ;; collect the bins
   binGrid = [1E,2E]
   binsizes = [0.3E,0.3E]
   binNames = ['In '+molecule,'Out '+molecule]
endif else begin
   if keyword_set(matchgrid) then begin
      startIndmatch = 12
      EndIndmatch = 478
      binGrid = lamgrid[startIndmatch:EndIndmatch]
      Nwavbins = n_elements(lamgrid[startIndmatch:EndIndmatch])
      binsizes = Dlam
   endif
   binGrid = (EndWav - StartWav) * findgen(Nwavbins)/float(Nwavbins) + $
             StartWav
   binsizes = fltarr(Nwavbins) + (EndWav - StartWav)/float(Nwavbins)
   
   if keyword_set(nyquist) then begin
      assert,Nwavbins,'>',1,"Warning, must be more than 1 bin for Nyquist sample"
      ;; find the bins in between with same bin width
      InbetweenBins = BinGrid[lindgen(Nwavbins-1l)] + binsizes[lindgen(Nwavbins-1l)]/2E
      InbetweenSizes = binsizes[lindgen(Nwavbins-1l)]
      ;; increase the number of bins
      Nwavbins = Nwavbins * 2l - 1l
      binsizes = [binsizes,InbetweenSizes]
      combinedGrid = [binGrid,InbetweenBins]
      binGrid = combinedGrid[sort(combinedGrid)]
   endif
endelse

binfl = fltarr(Nwavbins,nfile) ;; binned flux ratio
binflE = fltarr(Nwavbins,nfile)
binind = fltarr(Nwavbins,Nap,nfile) ;; individual binned fluxes
binindE = fltarr(Nwavbins,Nap,nfile)

if keyword_set(molecbin) then begin
   for i=0l,nfile-1l do begin
;      binfl[0,i] = 
;      binflE[*,i] = robust_sigma(Divspec
      if total(finite(Divspec[molecInd,0,i])) GT 5 then begin
         resistant_mean,Divspec[molecInd,0,i],10,mean1
         binfl[0,i] = mean1
         binflE[0,i] = robust_sigma(Divspec[molecInd,0,i])
         resistant_mean,Divspec[outmolec,0,i],10,mean2
         binfl[1,i] = mean2
         binflE[1,i] = robust_sigma(Divspec[molecInd,0,i])
      endif
      for k=0l,Nap-1l do begin
         if total(finite(flgrid[molecInd,k,i])) GT 5 then begin
            resistant_mean,flgrid[molecInd,k,i],10,mean1
            binind[0,k,i] = mean1
            binindE[0,k,i] = robust_sigma(flgrid[molecInd,k,i])
            resistant_mean,flgrid[outmolec,k,i],10,mean2
            binind[1,k,i] = mean2
            binindE[1,k,i] = robust_sigma(flgrid[outmolec,k,i])
         endif
      endfor
   endfor
endif else begin
   if keyword_set(matchgrid) then begin
      binfl[*,*] = Divspec[startIndmatch:EndIndmatch,0,*]
      binflE[*,*] = Divspec[startIndmatch:EndIndmatch,0,*]/SNR[startIndmatch:EndIndmatch,0,*]
      binind[*,*,*] = flgrid[startIndmatch:EndIndmatch,*,*]
      binindE[*,*,*] = Errgrid[startIndmatch:EndIndmatch,*,*]
   endif else begin
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
   endelse
endelse


if keyword_set(psplot) then begin
   device,/close
   cgPS2PDF,plotnm
   ;; return the plotting to normal
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

if not keyword_set(skipBJD) then begin
   ;; Convert the JD data to Barycentric coordinates

   BJD_TDB = UTC2BJD(utgrid,raDeg[0],decDeg[0],earthobs=obscode[0])
   utgrid = BJD_TDB

endif

; Save all data
save,flgrid,lamgrid,utgrid,bingrid,binfl,binflE,$
     ErrGrid,SNR,Divspec,DivspecE,backgrid,$
     Nwavbins,binsizes,binind,binindE,filen,$
     airmass,altitude,backgrid,header,$
     filename='data/specdata.sav'
end
