pro gather_mcmc_radii
;; collects the output of the MCMC results and puts it in an single
;; radius vs. wavelength file

  cd,c=currentd
  fileopt = file_search(currentd+'/data/mcmc/param_unc/param_unc*um.txt')
  nfile = n_elements(fileopt)

  
  for i=0l,nfile-1l do begin
     readcol,fileopt[i],paramNames,lmfit,limfitErr,mcmcfit,mcmcfitErrP,mcmcfitErrM,$
             format='(A,F,F,F,F,F)',skipline=1,/silent
     if i EQ 0l then begin
        nparam = n_elements(paramNames)
        ParamArray = fltarr(nfile,nparam)
        ParamArrayErr = fltarr(nfile,nparam)
     endif
     ParamArray[i,*] = mcmcfit
     ParamArrayErr[i,*] = (mcmcfitErrP + mcmcfitErrM)/2E
  endfor

  ;; Get the wavelengths from the Rp/R* file
  radfile='radius_vs_wavelength/radius_vs_wavl.txt'
  readcol,radfile,wavl,wavlsize,rad,rade,skipline=1,format='(F,F,F)'
  assert,n_elements(wavl),'=',nfile,'Mismatch between rad vs wavelength wavelengths and input wavelengths.'

  ;; Save the radius file
  forprint,wavl,wavlsize,ParamArray[*,0],ParamArrayErr[*,0],$
           comment='#Wavelength(um) Bin size(um) Rp/R*   Rp/R* Error',$
           textout='radius_vs_wavelength/mcmc_rad_vs_wavl.txt'

end

