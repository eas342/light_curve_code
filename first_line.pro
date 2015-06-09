function first_line,filen
if file_exists(filen) then begin
   
   st = ''
   
   openr,1,filen
   readf,1,st
   close,1
   return,st
endif else begin
   message,'File does not exist',/cont
   return,''
endelse

end
