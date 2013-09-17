pro compile_phot,dec23=dec23,dec29=dec29
;; Compiles the MORIS photometry data in the same way as spectra for
;; use by other scripts
;; dec23 -- look at the dec23 data set (default is jan04)
;; dec29 -- look at the dec29 data set (default is jan04)

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


Nwavbins = 1 ;; photometry!
bingrid = [0.872] ;; Z band photometry from curve
binsizes = [0.076]

npoints = n_elements(bjd)
utgrid = bjd + 2400000.5D ;- 0.0055D

binfl = fltarr(Nwavbins,npoints) ;; binned flux ratio
binflE = fltarr(Nwavbins,npoints)
binfl[0,*] = best_flux
binflE[0,*] = 1E

airmass = fltarr(npoints) + 1.E
altitude = fltarr(npoints) + 90E

morisPhase = phase

wavname='z-prime'

; Save all data the same way as spec
save,bingrid,binfl,binflE,$
     Nwavbins,binsizes,$
     airmass,altitude,$
     utgrid,morisPhase,wavname,$
     filename='data/specdata.sav'
end
