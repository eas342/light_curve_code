pro compile_spec,extraction2=extraction2,sum=sum,nwavbins=nwavbins,$
                 dec23=dec23,dec29=dec29,nyquist=nyquist,extremeRange=extremeRange,$
                 maskwater=maskwater,custRange=custRange,widewatermask=widewatermask,$
                 cleanbyeye=cleanbyeye,specshift=specshift,sameshift=sameshift,$
                 starshift=starshift,$
                 custmask=custmask,molecbin=molecbin,trycurved=trycurved,$
                 matchgrid=matchgrid,readCurrent=readCurrent,skipBJD=skipBJD,$
                 masktelluric=masktelluric,showall=showall,irafnoise=irafnoise,$
                 longwavname=longwavname,trycorrect=trycorrect,removelinear=removelinear,$
                 backRatio=backRatio,alreadyDivided=alreadyDivided,saveshifts=saveshifts,$
                 wavWeights=wavWeights,indbin=indbin,wavpixel=wavpixel,flipstars=flipstars,$
                 spatialRatio=spatialRatio,pretendTransit=pretendTransit,$
                 specKey=specKey,normalize=normalize,inject=inject,$
                 quickread=quickread
;; Compiles the spectra into a few simple arrays to look at the spectrophotometry
;; extraction2 -- uses whatever spectra are in the data directory
;; sum -- uses the variance weighted (optimal) extraction by
;;            default, but with the sum keyword it will use the
;;            standard sum extraction from the IRAF reduction
;; spatialRatio - uses the polynomial fit to the spatial ratio of the
;;                two stars instead of the sum of each star
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
;; sameshift -- force both stars to have the same wavelength shifts applied
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
;; showall - show all spectral info
;; irafnoise -- use IRAF to calculate noise instead of my own method
;; longwavname -- use longer wavelength name (describe wavelength range)
;; trycorrect -- tries to correct the flux by the background ratio
;; removelinear -- remove linear trends in the spectra before binning
;; backRatio -- use the ratio of backgrounds as a time series instead
;;              of ratio of fluxes
;; alreadyDivided -- the spectra are already divided so ignore the
;;                   DIVISOR keyword
;; SaveShifts -- saves the shifts instead of trying to correct for them
;; wavWeights -- weight the wavelength bins?
;; indBin -- bin the wavelengths for each star, then divide
;; wavpixel -- use the pixel numbers instead of wavelength for the
;;             wavelength grid
;; flipstars -- flip the reference and target stars
;; pretendTransit - shift the time to simulate what out-of-transit
;;                  data looks like in transit
;; specKey - allows one to set the spectral key by hand
;; normalize - normalize all flux ratios to be 1
;; quickread - re-read the FITS spectral images. Otherwise it uses a
;;          saved version

;Nwavbins = 35 ;; number of wavelength bins
;Nwavbins = 9 ;; number of wavelength bins
;Nwavbins = 20 ;; number of wavelength bins

if n_elements(Nwavbins) EQ 0 then Nwavbins = 9

SigRejCrit = 3 ;; number of sigma to reject when binning

if n_elements(wavWeights) EQ 0 then wavWeights=1

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
   readcol,'file_lists/current_speclist.txt',filen,format='(A)',stringskip='#'
   ;; in the case that your in the "readcurrent" mode, it's
   ;; probably with the latest pipeline that uses the IRAF-calculated
   ;; noise for optimal extraction. Still, if the use specifies that
   ;; the irafnoise keyword should be zero, then it will be left as such
   if n_elements(irafnoise) EQ 0 then irafnoise=1
endif
nfile = n_elements(filen)

;; Get the detector info
readcol,'data/detector_info.txt',skipline=1,$
       descrip,detectdata,format='(A,F)'
Gain = detectdata[0]
ReadN = detectdata[1]
;Npix = 14E ;; 14 pixel aperture

;;Get the coordinate info for the stars on the slit
readcol,'data/object_coordinates.txt',skipline=1,$
        objname,raHr,decDeg,format='(A,D,D)'
raDeg = raHr / 24D * 360D

;; get the planet info
readcol,'transit_info/planet_info.txt',info,data,format='(A,D)',$
        skipline=1
planetdat = create_struct('null','')
for l=0l,n_elements(info)-1l do begin
   planetdat = create_struct(planetdat,info[l],data[l])
endfor

;;Get the observatory
readcol,'data/observatory_info.txt',obscode,obsnm,$
        format='(A,A)',skipline=1

;; Get the length of 
a = mrdfits(filen[0],0,header)
sizea = size(a)
Ngpts = sizea[1] ;; number of grid points

;; Make a grid for the wavelengths
if keyword_set(wavpixel) then begin
   lamgrid = lindgen(Ngpts)
endif else begin
   if string(fxpar(header,"BANDID6")) EQ "Wavelength" then begin
      lamgrid = a[*,0,5]
   endif else begin
      DLam = fxpar(header,"CD1_1")
      lamstart = fxpar(header,"CRVAL1")
      ;; For IRAF dispersion corrected files
      lamgrid = (DLam * findgen(Ngpts) + lamstart)/1E4 ;; in microns
   endelse
endelse

;; Make a grid for the flux, background and UT time
Nap = sizea[2] ;; number of apertures
flgrid = fltarr(Ngpts,Nap,nfile)
backgrid = fltarr(Ngpts,Nap,nfile)
errgrid = fltarr(Ngpts,Nap,nfile)
if keyword_set(spatialRatio) then begin
   ratioGrid = fltarr(Ngpts,1,nfile)
endif
utgrid = dblarr(nfile)
itimegrid = dblarr(nfile)
airmass = dblarr(nfile)
Altitude = dblarr(nfile,Nap) ;; Altitude of each star
focus = fltarr(nfile)

case 1 of 
   n_elements(specKey) NE 0: specKey=specKey
   keyword_set(sum): SpecKey = 0
   keyword_set(BackRatio): SpecKey = 2
   else: SpecKey=1
endcase

;; Set up the aperture keys
ApKey = round([planetdat.TargStarNum[0],planetdat.RefStarNum[0]])

;; Since the H2RG detector, the array readout is flipped, so after
;; this time, I'm going to flip which star counts as #1 and #2

firstUTpt = double(fxpar(header,'MJD_OBS'))
;; The old data uses Julian date, but new uses MJD, so we have to add
;; the JD − 2400000.5
if firstUTpt LT 2400000D then firstUTpt = firstUTpt + 2400000.5D

ArrayChangeDate = date_conv('2014-08-01T00:00.00','JD')

if firstUTpt GT arrayChangeDate then begin
   apkey = reverse(apkey)
endif 
;; also have the ability to flip apertures
if keyword_set(flipstars) then begin
   apkey = reverse(apkey)
   print,"Flipping star order MANUALLY!"
endif

restore,'data/used_date.sav'
specDataCompressed = 'data/save_specdata/'+specfileListNamePrefix+'_specdata.sav'
if keyword_set(quickread) and file_exists(specDataCompressed) then begin
   prevNwavbins= nwavbins
   restore,specDataCompressed
   undefine,nwavbins
   nwavbins=prevNwavbins
   maxfile=-1
   prevFileUsed=1
endif else begin
   maxfile=nfile
   prevFileUsed=0
endelse

for i=0l,maxfile-1l do begin
   ;; Read all files into the grid
   a2 = mrdfits(filen[i],0,header2,/silent)
   ;; Check that the number of points is corrrect
   sizea2 = size(a2)

   if keyword_set(showall) then begin
      allcolors=mycol(['white','blue','green','yellow','white'])
      nalllines = 4
      allstyles=[0,2,3,4]

      for j=0l,Nap-1l do begin
         for k=0l,nalllines-1l do begin
            if k EQ 0 and j EQ 0 then begin
               plot,lamgrid,a2[*,j,k],color=allcolors[k],linestyle=allstyles[j],$
                    xtitle='Wavelength (um)',ytitle='Flux (DN)',$
                    xrange=[0.6,2.5],xstyle=1,$
                    yrange=[0,1E5]
            endif else oplot,lamgrid,a2[*,j,k],color=allcolors[k],linestyle=allstyles[j]
         endfor
      endfor
      legend,[fxpar(header2,'BANDID1'),$
              fxpar(header2,'BANDID2'),$
              fxpar(header2,'BANDID3'),$
              fxpar(header2,'BANDID4'),$
              'Aperture 2'],color=allcolors,linestyle=[0,0,0,0,2]
   endif

   ;; Get the # non-destructive reads IMPORTANT b/c flux must be divided by NDR
   NDR = double(fxpar(header2,'NDR'))
   if keyword_set(irafnoise) or keyword_set(alreadyDivided) then begin
      divisor=1.0E
   endif else divisor = NDR


;   divisor = divisor / 1.5E ;; b/c of fowler sampling w/ n_max reads
;   divisor = 1.0E

   for j=0,Nap-1 do begin
      flgrid[*,j,i] = a2[*,ApKey[j],SpecKey] * Gain / divisor
      backgrid[*,j,i] = a2[*,ApKey[j],2] * Gain / divisor;; multiply by gain divide by divisor
      errgrid[*,j,i] = a2[*,ApKey[j],3] * Gain / divisor
   endfor
   if keyword_set(spatialRatio) then ratiogrid[*,0,i] = a2[*,0,6]

;; CORRECTION FACTOR - the errors must be scaled for a variety of
;;                     reasons:
;; 1) The actual # of photons is smaller because SpeX gives the sum of
;; non-destructive reads, not the average of reads
;; 2) There are subtleties to the Fowler read process where the reads
;; are correlated with each other
;; see Reasearch Notes XII, page 52 for details
   itimeGrid[i] = double(fxpar(header2,'ITIME')) ;; (sec)
   Teff = itimeGrid[i]
   rtime = double(fxpar(header2,'TABLE_MS'))/1000E;; read time (sec)
   Tint = Teff + NDR * rtime
   aparray = double(strsplit(fxpar(header2,'APNUM1'),' ',/extract)) ;; array
   npix = aparray[3] - aparray[2] ;; number of pixels extracted
   nmax = (Tint)/(2E * rtime) ;; max read possible
   eta = NDR / nmax ;;read time duty cycle
   ReadN = 12E * sqrt(32E/NDR)

;   CorFac = sqrt(NDR * Lparam) * sqrt(flgrid[*,*,i] + backgrid[*,*,i] + readN^2 * npix/Lparam)/$
;            (alpha * (1E - eta/2E) * $
;             sqrt(flgrid[*,*,i] + backgrid[*,*,i] + readN^2 * npix/NDR))
;   errgrid[*,*,i] = CorFac * errgrid[*,*,i]
   ;; Terminology straight from Garnett & Forrest 1993
   FluxGF = flGrid[*,*,i] / Teff
   if keyword_set(backratio) then begin
      ;; If we're going to call the flux the same as the
      ;; background flux, then the actual background flux will be zero
      BackGF = 0E * backGrid[*,*,i]
   endif else BackGF = backGrid[*,*,i] / Teff
   Signal = FluxGF * Teff
   NoisePhoton = sqrt((FluxGF + BackGF) * Tint * (1E - 2E * eta/3E + 1E/(6E * eta * nmax^2)))

   NoiseRead = sqrt(ReadN^2 * npix)
   TotNoise = sqrt(NoisePhoton^2 + NoiseRead^2)
   SignalNaive = FluxGF * NDR * Teff
   NoiseNaive = sqrt((FluxGF + BackGF) * Teff * NDR + npix * readN^2) ;; doesn't know about NDRs
;   CorFac = (TotNoise/Signal)/(NoiseNaive/SignalNaive)
;   flgrid[*,*,i] = fluxGF * Tint
;   backgrid[*,*,i] = backGF * Tint
   NoiseTry = sqrt((FluxGF + BackGF)* Tint + npix * ReadN^2 )

   if not keyword_set(irafnoise) then begin
      errGrid[*,*,i] = TotNoise
   endif

   utgrid[i] = double(fxpar(header2,'MJD_OBS'))
   utgrid[i] = utgrid[i] + itimeGrid[i]/(3600D * 24D)

   airmass[i] = double(fxpar(header2,'AIRMASS',count=airmassFinds))
   if airmassFinds EQ 0 then airmass[i] = double(fxpar(header2,'TCS_AM',count=airmassFinds))
   focus[i] = float(fxpar(header2,'TCS_FOC',count=focusFinds))
   ;; find the differential airmass between the two stars
   utarray = dblarr(Nap) + utgrid[i]
   eq2hor,raDeg,decDeg,utarray,alt,az,obsname=obscode[0]
   for j=0,Nap-1 do begin
      altitude[i,j] = alt[j]
   endfor
endfor

;; The old data uses Julian date, but new uses MJD, so we have to add
;; the JD − 2400000.5

if utgrid[0] LT 2400000D then utgrid = utgrid + 2400000.5D

;; For the out-of-transit KIC 1255, we can move the time to see what
;; out-of-transit data looks like when fit with transit models
if keyword_set(pretendTransit) then utgrid = utgrid - 0.30D

;; Reset all zeros and negative flux values
badp = where(flgrid LE 0)
if badp NE [-1] then flgrid[badp] = !values.f_nan

badback = where(backgrid LE 0)
if badback NE [-1] then backgrid[badback] = !values.f_nan

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
;   telluricW = [[1.1,1.15],[1.3,1.48],[1.78,1.95],[1.97,2.05]]
   telluricW = [[1.1,1.15],[1.3,1.48],[1.75,2.1]]
   sizeTell = size(telluricW)
   Ntelluric = sizeTell[2]
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

;; Factor for Fowler sampling
;FowFac = 1.0E ;; from Garnett & Forrest http://adsabs.harvard.edu/abs/1993SPIE.1946..395G
;; In the terminology of Garnett & Forrest, eta = 1
;ErrGrid = nansqrt( FowFac * flgrid + FowFac * backgrid + readnarr^2 * npix)


;; Star Shift, the default is to move the reference star by -1 pixel
if n_elements(starshift) EQ 0 then begin
   ;; the readCurrent keyword is designed to read newer data so the default star
   ;; shift is 0 for the new data
   if keyword_set(readCurrent) then starshift=0 else begin
      starshift = -1
   endelse
endif
xygrid = fltarr(Ngpts,nfile)
xygrid[*,*] = flgrid[*,1,*] ;; adjust the reference star
shiftedGrid = shift_interp(xygrid,starshift)
flgrid[*,1,*] = shiftedGrid

;; Shift arrays
if keyword_set(specshift) or keyword_set(saveShifts) then begin
   ;; Align the stars within their bins

   if keyword_set(trystraight) then begin
      ;; since the straightend spec has some weird properties above
      ;; 2.47um because of shifting spectra, mask those out
      badp = where(lamgrid GT 2.47)
      flgrid[badp,*,*] = !values.f_nan
   endif

   specShiftArr = fltarr(Nap,nfile)
   masterspec = median(flgrid[*,0,*],dimension=3)
   for i=0l,Nap-1l do begin
      xyspec = fltarr(Ngpts,nfile)
      xyspec[*,*] = flgrid[*,i,*]

      ShiftedGrid = find_shifts(xyspec,/cutEnds,masterspec=masterspec)
;      ShiftedGrid = find_shifts(xyspec,/cutEnds)
      restore,'data/wavelength_shifts/temp_shift_list.sav'
      specshiftArr[i,*] = shiftArr

;      if keyword_set(specshift) then flgrid[*,i,*] = shiftedgrid
   endfor
   ;; Shift the two stars by their average, with an average offset
   if keyword_set(specshift) then begin
      avgShifts = total(specshiftarr,1)/2E
      medDiff = median(specshiftarr[1,*] - specshiftarr[0,*])

      bestshifts = fltarr(Nap,nfile) ;; well maybe not the best after all!
      bestshifts[0,*] = avgShifts
      bestshifts[1,*] = avgShifts + medDiff

      for i=0l,Nap-1l do begin
         for j=0l,nfile-1l do begin
            if keyword_set(sameshift) then begin
               flgrid[*,i,j] = shift_interp(flgrid[*,i,j],bestshifts[i,j])
            endif else begin
               flgrid[*,i,j] = shift_interp(flgrid[*,i,j],specshiftarr[i,j])
               if keyword_set(spatialRatio) and i EQ 0 then begin
                  ratiogrid[*,i,j] = shift_interp(ratiogrid[*,i,j],specshiftarr[i,j])
               endif
            endelse
         endfor
      endfor
   endif

   ;; Save the shift arr
   if keyword_set(saveShifts) then begin
      ;; also get the speclist name
      restore,'data/used_date.sav'
      save,specshiftArr,filename='data/shift_data/shift_'+specfileListNamePrefix+'.sav'
   endif

endif

if keyword_set(trycorrect) then begin
   flgrid[*,0,*] = median(flgrid[*,0,*]) * 3E + flgrid[*,0,*]
endif

if keyword_set(spatialRatio) then begin
   Divspec = ratioGrid
endif else begin
;; Divide the two spectra
   Divspec = flgrid[*,0,*] / flgrid[*,1,*]
endelse

if n_elements(inject) GT 0 then begin
   inject_transit,float(inject),divspec=divspec,utgrid=utgrid
endif

if keyword_set(removelinear) then begin
   remove_linear,utgrid,divSpec,lamgrid
endif
if keyword_set(normalize) then begin
   medianNormArrDiv = rebin(median(divspec,dimension=3),Ngpts,1,nfile)
   divspec = divspec/medianNormArrDiv
endif

fracE = nansqrt((ErrGrid[*,0,*]/flgrid[*,0,*])^2 + $
             (ErrGrid[*,1,*]/flgrid[*,1,*])^2 )

DivspecE = fracE * Divspec
SNR = Divspec / DivspecE

;; Divide the two backgrounds
backdiv = backgrid[*,0,*] / backgrid[*,1,*]

if keyword_set(trycorrect) then begin
;   C0 = 1.15E
;   C1 = -0.2E
   Divspec = (median(flgrid[*,0,*]) * 0E + flgrid[*,0,*]) / (flgrid[*,1,*])

   DivspecE = fracE * Divspec
endif


;; Do the wavelength binning for the divided spec
if keyword_set(extremeRange) then begin
   StartWav = lamgrid[0]
   EndWav = lamgrid[Ngpts-1]
endif else begin
   if n_elements(custRange) EQ 0 then begin
      readcol,'data/wavelength_ranges.txt',StartWavArr,EndWavArr,$
              /silent,format='(F,F)',skipline=1
      StartWav = StartWavArr[0]
      EndWav = EndWavArr[0]
   endif else begin
      StartWav = custRange[0]
      EndWav = custRange[1]
   endelse
endelse
;; These are the starts of the bins (not the middles)
if keyword_set(molecbin) then begin
   Nwavbins = 2 ;; in the molecule and outside the molecule
   molecule='H2O'
   readcol,'molec_inputs/H2O_bin_locations.txt',molecStarts,molecEnds,$
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
   wavName = ['In-'+molecule,'Out-'+molecule]
endif else begin
   if keyword_set(matchgrid) then begin
      tabinv,lamgrid,startWav,startIndmatch
      startIndmatch = round(startIndMatch)
      tabinv,lamgrid,endWav,EndIndmatch
      binGrid = lamgrid[startIndmatch:EndIndmatch]
      Nwavbins = n_elements(lamgrid[startIndmatch:EndIndmatch])
      binsizes = lamgrid[1]-lamgrid[0]
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
binback = fltarr(Nwavbins,Nap,nfile) ;; individual binned background flux
binbackE = fltarr(Nwavbins,Nap,nfile) ;; individual binned background flux error


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
         y = avg_series(lamgrid,Divspec[*,0,i],SNR[*,0,i],binGrid,binsizes,weighted=wavWeights,$
                        oreject=sigRejCrit,eArr=yerr,/silent,errIn=divSpecE[*,0,i])
         binfl[*,i] = y
         binflE[*,i] = yerr

         for k=0,Nap-1 do begin
            y2 = avg_series(lamgrid,flgrid[*,k,i],flgrid[*,k,i]/Errgrid[*,k,i],binGrid,$
                            binsizes,weighted=wavWeights,$
                            oreject=sigRejCrit,eArr=yerr2,/silent,errIn=Errgrid[*,k,i])
            binind[*,k,i] = y2
            binindE[*,k,i] = yerr2
            y3 = avg_series(lamgrid,backgrid[*,k,i],fltarr(Ngpts)+1E,binGrid,$
                            binsizes,weighted=wavWeights,$
                            oreject=sigRejCrit,/silent,stdevArr=binbackStDev)

            binback[*,k,i] = y3
            binbackE[*,k,i] = binbackStDev
         endfor
      endfor
   endelse
endelse

if keyword_set(indbin) then begin
   ;; Bin each stars first, then divide the two stars instead of
   ;; dividing the two stars, then binning
   binfl = fltarr(Nwavbins,nfile) ;; binned flux ratio
   binflE = fltarr(Nwavbins,nfile)
   for k=0l,Nwavbins-1l do begin
      goodp = where(binind[k,1,*] NE 0)
      binfl[k,goodp] = binind[k,0,goodp]/binind[k,1,goodp]
      binflE[k,goodp] = binfl[k,goodp] * sqrt($
                        (binindE[k,0,goodp]/binind[k,0,goodp])^2 +$
                        (binindE[k,1,goodp]/binind[k,1,goodp])^2)
   endfor
endif


if keyword_set(psplot) then begin
   device,/close
   cgPS2PDF,plotnm
   ;; return the plotting to normal
   device,decomposed=0
   set_plot,'x'
   !p.font=-1
endif

if not keyword_set(skipBJD) and not prevFileUsed then begin
   ;; Convert the JD data to Barycentric coordinates
   BJD_TDB = UTC2BJD(utgrid,raDeg[0],decDeg[0],earthobs=obscode[0])
   utgrid = BJD_TDB

endif

;; Describe the wavlengths
bingridmiddle = bingrid + binsizes/2E
if not keyword_set(molecbin) then begin
   wavname = strarr(Nwavbins)
   if keyword_set(wavpixel) then wavFormat='(F04.0)' else wavFormat='(F4.2)'
   for k=0l,Nwavbins-1l do begin
      if keyword_set(longwavname) then begin
         wavname[k] = string(bingrid[k],format=wavFormat)+'-'+$
                      string(bingrid[k]+binsizes[k],format=wavFormat)
      endif else wavname[k] = string(bingridmiddle[k],format=wavFormat)
   endfor
endif

; Save all data
save,flgrid,lamgrid,bingrid,binfl,binflE,backdiv,$
     utgrid,itimegrid,wavname,$
     ErrGrid,SNR,Divspec,DivspecE,backgrid,$
     Nwavbins,binsizes,binind,binindE,filen,$
     airmass,altitude,header,apkey,binback,binbackE,$
     focus,$
     filename='data/specdata.sav'

file_copy,'data/specdata.sav',specDataCompressed,/overwrite

end
