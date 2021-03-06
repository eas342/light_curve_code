pro change_kernels,option
;; Switches the type of covariance kernel used in MCMC fitting

if n_elements(option) EQ 0 then begin
   filen = choose_file(searchDir='data/kernels',filetype='.txt')
endif else filen=option

cd,c=currentd
kernelfileN = currentd+'/cov_kernel.pro'

spawn,'cp '+filen+' '+kernelfileN

resolve_routine,'cov_kernel',/is_function

end
