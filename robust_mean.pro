function robust_mean,array,oreject=oreject,err=err
;; Does outlier rejection and calculates the mean of an array
;; Also calculates the error in the mean


sig = robust_sigma(array)

med = median(array)

if n_elements(oreject) EQ 0 then goodp = where(finite(array)) else begin
   goodp = where(finite(array) and $
                 array LT med + double(oreject) * sig and $
                 array GT med - double(oreject) * sig,ngood)
endelse

if goodp EQ [-1] then return,!values.d_nan else begin
   meanV = mean(array[goodp])

   err = sig / sqrt(double(ngood) -1D)

   return, meanV
endelse

end
