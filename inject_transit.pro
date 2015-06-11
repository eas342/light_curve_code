pro inject_transit,strength,divspec=divspec
;; Artificially injects a transit into the control night

if n_elements(strength) EQ 0 then strength=1E
;; Get the data

if n_elements(divspec) EQ 0 then begin
   restore,'data/specdata.sav' ;; get the spectro-photometry data
   DoSave=1
endif else doSave=0
;restore,'data/timedata.sav' ;; get the orbital phase

tplot = find_phase(utgrid)

sz = size(divspec)

ntime = n_elements(tplot)

;; Change the individual divspec image

KepTransit = kepler_func(tplot,strength)
modImg = rebin(reform(KepTransit,[1,1,ntime]),sz[1],sz[2],sz[3])
divspec = divspec * modImg

;; Change the time series image

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
