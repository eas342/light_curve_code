pro compile_phot,dec23=dec23,dec29=dec29,readC=readC,removelinear=removelinear,$
                 choosefile=choosefile,thphot=thphot,inject=inject,$
                 pretendTransit=pretendTransit
;; Compiles the MORIS photometry data in the same way as spectra for
;; use by other scripts
;; dec23 -- look at the dec23 data set (default is jan04)
;; dec29 -- look at the dec29 data set (default is jan04)
;; readC - use the photometry for kic1255
;; choosefile -- choose a specific file from the CoRoT-1 photometry
;;               folder
;; thphot - use the photometry that Terry gave me
;; inject - for the control night, inject a transit to see if it is recovered
;; pretendTransit - shifts the time series to simulate the fitting
;;                  process on a transit 

case 1 of 
   keyword_set(choosefile): begin
      custfile = choose_file(searchDir='../moris_data/reduced_lightc/',filetype='.sav')
      restore,custfile
   end
   keyword_set(readC): begin
      restore,'../../../Documents/kic1255/extracted_photometry/zhao_photometry/kic1255_npoly2.idl'
      restore,'data/used_date.sav'
      case useDate of 
         '2013aug13': datastruct = data0813
         '2013aug15': datastruct = data0815
         '2013aug17': datastruct = data0818 ;; mistyped in the photometry save file
         '2013sep03': datastruct = data0903
         '2014sep02': begin
            restore,'../moris_data/reduced_lightc/kic1255_tlh_UT2014.09.02_data.sav'
         end
         '2014sep03': begin
            restore,'../moris_data/reduced_lightc/kic1255_tlh_UT2014.09.03_data.sav'
         end
         '2014sep04': begin
            restore,'../moris_data/reduced_lightc/kic1255_tlh_UT2014.09.04_data.sav'
         end
         '2014aug14': begin
            restore,'../moris_data/reduced_lightc/kic1255_tlh_UT2014.08.14_data.sav'
         end
         '2014aug18': begin
            restore,'../moris_data/reduced_lightc/kic1255_tlh_UT2014.08.18_data.sav'
         end
      endcase
      bjd = datastruct.bjd_tdb
      best_flux = datastruct.flux
   end
   keyword_Set(thphot): begin
      restore,'../moris_data/reduced_lightc/kic1255_tlh_UT2014.09.04_data.sav'
   end
   else: begin
      case 1 of
         keyword_set(dec23): begin
            restore,'../moris_data/reduced_lightc/corot-1_UT2011Dec23_a0729_final.idl'
         end
         keyword_set(dec29): begin
            restore,'../moris_data/reduced_lightc/UT2011Dec29_corot-1_a0729_final.idl'
         end
         else: begin
            restore,'../moris_data/reduced_lightc/UT2012Jan04_corot-1_a0729_final.idl'
         end
      endcase
   end
endcase


Nwavbins = 1 ;; photometry!
if keyword_set(readC) or keyword_set(thphot) then begin
   bingrid = [0.626] ;; R band photometry from curve
   binsizes = [0.070]
   wavname='r-prime'
   utgrid = bjd
endif else begin
   bingrid = [0.826] ;; Z band photometry from curve
   binsizes = [0.076]
   wavname='z-prime'
   utgrid = bjd + 2400000.5D    ;- 0.0055D
endelse


if keyword_set(pretendTransit) then utgrid = utgrid - 0.3

npoints = n_elements(bjd)

binfl = fltarr(Nwavbins,npoints) ;; binned flux ratio
binflE = fltarr(Nwavbins,npoints)
binfl[0,*] = best_flux
if keyword_set(readC) then begin
   binflE[0,*] = datastruct.err
endif else binflE[0,*] = 1E

airmass = fltarr(npoints) + 1.E
altitude = fltarr(npoints) + 90E

if n_elements(phase) NE 0 then morisPhase = phase

if n_elements(inject) GT 0 then begin
   inject_transit,float(inject),divspec=binfl,utgrid=utgrid,$
                  /photMode
endif

if keyword_set(removelinear) then begin
   remove_linear,utgrid,binfl,bingrid
endif

; Save all data the same way as spec
save,bingrid,binfl,binflE,$
     Nwavbins,binsizes,$
     airmass,altitude,$
     utgrid,morisPhase,wavname,$
     filename='data/specdata.sav'
end
