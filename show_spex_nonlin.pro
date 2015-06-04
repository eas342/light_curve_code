pro show_spex_nonlin
;; Plots the non-linearity of the Aladin InSb SpeX detector

restore,filepath('lc_coeff.sav', ROOT=file_dirname(file_which('SpeX.dat'),/MARK))
;; Gets lc_coeff

imgSizeX = 1024
imgSizeY = 1024

nstep = 10
fullL = 8500E
Gain=12.1E
countArr = findgen(nstep)/float(nstep-1) * fullL
tempimg = fltarr(imgSizeX,imgSizeY)


medCor = fltarr(nstep)

for i=0l,nstep-1l do begin
   Corr = mc_imgpoly(tempimg + countArr[i],lc_coeff)/lc_coeff[*,*,0]
   medCor[i] = median(Corr)
endfor
dat = struct_arrays(create_struct('INPUT',countArr,$
                                  'OUTPUT',medCor * countArr,$
                                  'NONLIN',medCor))
gparam = create_struct('TITLES',['Input Energy (DN)','Detected Energy (DN)',''],$
                      'LEGTITLE','','SLABEl',['SpeX Aladin Detector','Ideal'])
dat2 = dat                       
ev_oplot,dat,countArr,countArr,gparam=gparam
if keyword_set(overall) then begin
   genplot,dat,gparam=gparam
endif else begin
   gparam=create_struct('TITLES',['Input Energy (DN)','Non-linearity correction',''],$
                        'PKEYS',['INPUT','NONLIN'])
   genplot,dat2,gparam=gparam
endelse

end
