pro collect_spec,nwavbins=nwavbins
  filen = '../terry_reduction/CoRot1b_dsratio_detrend_v2.fits'
  utfilen = '../terry_reduction/utgrid_dec23.sav'
  ind1filen = '../terry_reduction/CoRot1b_ds_v2.fits'
  ind2filen = '../terry_reduction/CoRot1b_ds_v2.fits'

  restore,utfilen

  if n_elements(nwavbins) EQ 0 then nwavbins=9
  Nap =2
  wavweights = 0
  sigRejCrit = 5E

  data = mrdfits(filen,0,header)
  nfile = n_elements(data[0,*]) ;number of files


  startwav=0 ;; starting wavelength
  nwavs = n_elements(data[*,0])
  endWav = nwavs-1l
  lamgrid = lindgen(nwavs) + startwav


  flgrid = fltarr(nwavs,2,nfile)
  star1 = mrdfits(ind1filen,0,header)
  flgrid[*,0,*] = star1[*,*]
  star2 = mrdfits(ind2filen,0,header)
  flgrid[*,1,*] = star2[*,*]

  
  binGrid = (EndWav - StartWav) * findgen(Nwavbins)/float(Nwavbins) + $
            StartWav
  binsizes = fltarr(Nwavbins) + (EndWav - StartWav)/float(Nwavbins)
   

  binfl = fltarr(Nwavbins,nfile) ;; binned flux ratio
  binflE = fltarr(Nwavbins,nfile)
  binind = fltarr(Nwavbins,Nap,nfile) ;; individual binned fluxes
  binindE = fltarr(Nwavbins,Nap,nfile)


  ;; Make a divspec for reading in. Make the SNR uniform
  divspec = fltarr(Nwavs,1,nfile)
  divspec[*,0,*] = data[*,*]
  divspecE = fltarr(Nwavs,1,nfile) + 0.02E
  SNR = fltarr(Nwavs,1,nfile) + 50E

  for i=0,nfile-1 do begin
     y = avg_series(lamgrid,Divspec[*,0,i],SNR[*,0,i],binGrid,binsizes,weighted=wavWeights,$
                    oreject=sigRejCrit,eArr=yerr,/silent,errIn=divSpecE[*,0,i])
     binfl[*,i] = y
     binflE[*,i] = yerr
     
;     for k=0,Nap-1 do begin
;        y2 = avg_series(lamgrid,flgrid[*,k,i],flgrid[*,k,i]/Errgrid[*,k,i],binGrid,$
;                        binsizes,weighted=wavWeights,$
;                        oreject=sigRejCrit,eArr=yerr2,/silent,errIn=Errgrid[*,k,i])
;        binind[*,k,i] = y2
;        binindE[*,k,i] = yerr2
;     endfor
  endfor
;; Describe the wavlengths
  bingridmiddle = bingrid + binsizes/2E
  if not keyword_set(molecbin) then begin
     wavname = strarr(Nwavbins)
     for k=0l,Nwavbins-1l do begin
        if keyword_set(longwavname) then begin
           wavname[k] = string(bingrid[k],format='(F4.2)')+'-'+$
                        string(bingrid[k]+binsizes[k],format='(F4.2)')
        endif else wavname[k] = string(bingridmiddle[k],format='(F5.0)') + ' px'
     endfor
  endif
  
; Save all data
  ;; make the airmass flat for now
  airmass = fltarr(nfile) + 1.0E

  save,lamgrid,bingrid,binfl,binflE,utgrid,$;backdiv,$
       ;itimegrid
       wavname,divspec,flgrid,$
;       ErrGrid,SNR,Divspec,DivspecE,backgrid,$
       Nwavbins,binsizes,binind,binindE,$;filen,$
       airmass,header,$;altitude,backgrid,header,apkey,$
       filename='data/specdata.sav'
  
  
end
