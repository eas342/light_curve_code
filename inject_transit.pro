pro inject_transit,strength,divspec=divspec,$
                   photMode=photMode,utgrid=utgrid
;; Artificially injects a transit into the control night
;; divspec is the divided spectra from compile_spec
;; pretendTransit - passes the pretend transit on which is useful for
;;                  photometry mode
;; photMode - treats divspec as a binned flux for photometry
;;            (different size and dimension of array)
;; utgrid - if input, the JDB_UTC which needs to be converted to flux

if n_elements(strength) EQ 0 then strength=1E
;; Get the data

if n_elements(divspec) EQ 0 then begin
   restore,'data/specdata.sav' ;; get the spectro-photometry data
   DoSave=1
endif else doSave=0
;restore,'data/timedata.sav' ;; get the orbital phase

tplot = find_phase(utgrid,pretendTransit=pretendTransit)

sz = size(divspec)
ntime = n_elements(tplot)

;; Change the individual divspec image

KepTransit = kepler_func(tplot,strength)

if keyword_set(photMode) then begin
   modImg = reform(KepTransit,[1,ntime])
   divspec = divspec * modImg
endif else begin
   modImg = rebin(reform(KepTransit,[1,1,ntime]),sz[1],sz[2],sz[3])
   divspec = divspec * modImg
endelse



if DoSave then begin
; Save all data
   save,flgrid,lamgrid,bingrid,binfl,binflE,backdiv,$
        utgrid,itimegrid,wavname,$
        ErrGrid,SNR,Divspec,DivspecE,backgrid,$
        Nwavbins,binsizes,binind,binindE,filen,$
        airmass,altitude,backgrid,header,apkey,$
        focus,$
        filename='data/specdata.sav'
endif

end
