pro compile_multi_night

;; get the list of speclists
  readcol,'file_lists/multi_night.txt',listFile,format='(A)'

for i=0l,3l-1l do begin
   spawn,'cp '+listFile[i]+' file_lists/current_speclist.txt'

   ;; Get the spectral data, removing linear trends in the time series
   compile_spec,/readC,/removelinear 
   restore,'data/specdata.sav'

   if i EQ 0l then begin
      binflNew = binfl
      binflENew = binflE
      utgridNew = utgrid
   endif else begin
      binflNew = transpose([transpose(binflNew),transpose(binfl)])
      binflENew = transpose([transpose(binflENew),transpose(binflE)])
      utgridNew = [utgridNew,utgrid]
   endelse
   
endfor

binfl = binflNew ;; binned flux ratio
binflE = binflENew ;; standard error
utgrid = utGridNew

; Save all data
save,flgrid,lamgrid,bingrid,binfl,binflE,backdiv,$
     utgrid,itimegrid,wavname,$
     ErrGrid,SNR,Divspec,DivspecE,backgrid,$
     Nwavbins,binsizes,binind,binindE,filen,$
     airmass,altitude,backgrid,header,$
     filename='data/specdata.sav'
end
