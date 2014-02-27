pro compile_both,dec23=dec23,dec29=dec29,readC=readC,$
                 custwavbins=custwavbins,masktelluric=masktelluric,$
                 removelinear=removelinear,specshift=specshift
;; Compiles the both MORIS photometry & SpeX data so they can be
;; plotted & fit together
;; dec23 -- look at the dec23 data set (default is jan04)
;; dec29 -- look at the dec29 data set (default is jan04)
;; readC -- read the current speclist (from choose_speclist)
;; other commands are what's passed on to compile spec
;; specshift - passes this keyword onto compile_spec

;; Get the MORIS data
case 1 of
   keyword_set(dec23): compile_phot,/dec23,readC=readC,removelinear=removelinear
   keyword_set(dec29): compile_phot,/dec29,readC=readC,removelinear=removelinear
   else: compile_phot,readC=readC,removelinear=removelinear
endcase

restore,'data/specdata.sav'

nbands = Nwavbins
bandcenter = binGrid
bandwidth = binsizes
bandname = wavname
nPhot = n_elements(utgrid)

;; get the epoch info
;plot_tim_ser
;restore,'data/timedata.sav'

binPhot = binfl
binPhotE = binflE
timePhot = utgrid

case 1 of
   keyword_set(dec23): compile_spec,/dec23,readC=readC,nwavbins=custwavbins,$
                                    masktelluric=masktelluric,removelinear=removelinear,$
                                    specshift=specshift
   keyword_set(dec29): compile_spec,/dec29,readC=readC,nwavbins=custwavbins,$
                                    masktelluric=masktelluric,removelinear=removelinear
   else: compile_spec,readC=readC,nwavbins=custwavbins,$
                      masktelluric=masktelluric,removelinear=removelinear,$
                      specshift=specshift
endcase

;; Get the spectral data
restore,'data/specdata.sav'

nSpeX = n_elements(utgrid)
binflNew = fltarr(Nwavbins+nbands,nSpeX)
binflENew = fltarr(Nwavbins+nbands,nSpeX)
bingridNew = [bandcenter,bingrid]
binsizesNew = [bandwidth,binsizes]
wavnameNew = strarr(Nwavbins+nbands)

expDays = itimeGrid/(24D * 3600D)

;; Bin the photometry to fit in the same time bins as spectra
for i=0,nbands-1l do begin
   ;; Clean the photomery
;   ymodelIn = quadlc(tplot,planetdat.p,planetdat.b_impact,$
;                              0E,0E,planetdat.a_o_rstar)
   ybin = avg_series(timePhot,binPhot[i,*],binPhot[i,*]/binPhotE[i,*],$
                     utgrid - expDays,expDays,weighted=1,$
                     oreject=5E,/silent,errIn=binPhotE[i,*],stdevArr=stdevArr)
   binflNew[i,*] = ybin
   binflENew[i,*] = replicate(median(stdevArr),nSpeX)
   wavnameNew[i] = bandname[i]
endfor

for i=nbands,nbands+Nwavbins-1l do begin
   binflNew[i,*] = binfl[i-nbands,*]
   binflENew[i,*] = binfl[i-nbands,*]
   wavnameNew[i] = wavname[i-nbands]
endfor

binfl = binflNew ;; binned flux ratio
binflE = binflENew ;; standard error
bingrid = bingridNew
binsizes = binsizesNew
wavname = wavnameNew
Nwavbins = nbands+Nwavbins

; Save all data
save,flgrid,lamgrid,bingrid,binfl,binflE,backdiv,$
     utgrid,itimegrid,wavname,$
     ErrGrid,SNR,Divspec,DivspecE,backgrid,$
     Nwavbins,binsizes,binind,binindE,filen,$
     airmass,altitude,backgrid,header,$
     filename='data/specdata.sav'
end
