pro chisquare_compare,alt=alt,modchoose=modchoose,addbean=addbean,ironly=ironly,$
                      addalonso=addalonso
;; Calculates the chi-squared of spectral models with the data stored
;; in avg rp_rs
;; alt - use the alternate planet radius from Barge et al. 2008
;; addbean - add the Bean 2009 data
;; addalonso - add the Alonso et al. 2008 data
;; ironly - the data is only SpeX data - no photometry is included

readcol,'radius_vs_wavelength/avg_rp_rs.txt',wavl,wavlsize,rad,rade,$
        skipline=1,format='(F,F,F,F)'

if keyword_set(addbean) or keyword_set(addalonso) then begin
;; Add the Bean 2009 result
   CorWav = 0.65         ;; microns, approximately
   CorWid = 0.20         ;; microns, approx
   case 1 of
      keyword_set(alt): begin
         CorRad = 0.1388 ;;Rp/R*, Barge 2008
         CorErr = 0.0021
      end
      keyword_set(addalonso): begin
         CorRad = 0.1667 ;; Rp/R* Alonso et al. 2008
         CorErr = 0.0006
      end
      else: begin
         CorRad = 0.1433 ;;Rp/R*, Bean 2009
         CorErr = 0.0010
      end
   endcase
   wavl = [CorWav,wavl]
   wavlsize = [CorWid,wavlsize]
   rad = [CorRad,rad]
   rade = [CorErr,rade]
endif   

wavlwidth = wavlsize

if keyword_set(modchoose) then begin
   theofile = choose_file(searchDir='../models/transit_models/',filetype='.dat')
   readcol,theofile,theowav,theorad,$
           format='(F,F)'
endif else begin
;   readcol,'../models/fortney_g10mps_2500K_isothermal.csv',theowav,theorad,$
;           skipline=6,format='(F,F)'
   
   readcol,'../models/transit_models/lambda_R_P_iso_g10_2500.dat',theowav,theorad,$
           skipline=6,format='(F,F)'
   
endelse

;; Get the zprime filter curves
readcol,'../calculations/zprime_transmission/zprime_response.txt.csv',skipline=1,$
        zprimeWavl,zprimeResp
readcol,'../corot_data/filter_curve/filter_curve_CoRoT.txt',skipline=1,$
        CoRoTtransWav,CoRoTtrans
;; Change from nm to microns
CoRoTtransWav = CoRoTtransWav * 1E-3

ntheo=n_elements(theorad)
binModel1 = avg_series(theowav,theorad,fltarr(ntheo)+0.2E,wavl-wavlwidth/2E,wavlwidth,weighted=0)

if not keyword_set(ironly) then begin
;; Correct the photometry bins, starting with the CoRoT response
   binModel1[0] = photobin(theowav,theorad,CoRoTtransWav,CoRoTtrans)
   binModel1[1] = photobin(theowav,theorad,zprimeWavl,zprimeResp)
endif

binModel = binModel1
save,wavl,binModel,filename='data/binned_model.sav'

nwavs = n_elements(wavl)

start1 = [0.1D]
expr1 = 'model_evaluate(X,P[0])'
result = mpfitexpr(expr1,wavl,rad,rade,start1)
ymod1 = expression_eval(expr1,wavl,result)
dof1 = nwavs - n_elements(start1)
chisq1 = total(((rad - ymod1)/rade)^2)/dof1
print,'Literature Model Chisq/DOF= ',chisq1

start2 = [0.1D]
expr2 = 'P[0]'
resultFlat = mpfitexpr(expr2,wavl,rad,rade,start2)
ymod2 = replicate(resultFlat,nwavs)
dof2 = nwavs - n_elements(start2)
chisq2 = total(((rad - ymod2)/rade)^2)/dof2
print,'Flat Model ChisQ/DOF = ',chisq2

readcol,'../models/transit_models/transit_t2500g10_noTiO.dat',theowav2,theorad2,$
        skipline=6,format='(F,F)'
ntheo2 = n_elements(theorad2)

binModel2 = avg_series(theowav2,theorad2,fltarr(ntheo2)+0.2E,wavl-wavlwidth/2E,wavlwidth,weighted=0)
;; Correct the photometry bins, starting with the CoRoT response
binModel2[0] = photobin(theowav2,theorad2,CoRoTtransWav,CoRoTtrans)
binModel2[1] = photobin(theowav2,theorad2,zprimeWavl,zprimeResp)

binModel = binModel2
save,wavl,binModel,filename='data/binned_model.sav'

start3 = [0.1D]
expr3 = 'model_evaluate(X,P[0])'
result2 = mpfitexpr(expr3,wavl,rad,rade,start3)
ymod3 = expression_eval(expr3,wavl,result2)
dof3 = nwavs - n_elements(start3)
chisq3 = total(((rad - ymod3)/rade)^2)/dof3
print,'No TiO model Chisq/DOF = ',chisq3

plot,wavl,rad,/nodata,ystyle=16
oploterror,wavl,rad,wavlwidth/2E,rade,psym=3
oplot,wavl,ymod1,color=mycol('yellow')
oplot,wavl,ymod2,color=mycol('lblue')


;; Save all the data to file for use in plotting
nmod = 2

binnedWav = wavl
binnedValues = fltarr(nwavs,nmod)
binnedValues[*,0] = binModel1 * result[0]
binnedValues[*,1] = binModel2 * result2[0]

fullRes = create_struct('model1wav',theowav,'model1rad',theorad * result[0],$
                       'model2wav',theowav2,'model2rad',theorad2 * result2[0])
modName = ['Equilibrium Chemistry','TiO-Removed']

save,filename='data/binned_final.sav',$
     binnedValues,nmod,fullRes,binnedWav,$
     modName

end
