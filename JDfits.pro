function JDfits,filen,keyword
; takes a fits filename and keyword and returns the UT time associated
; with that keyword
hdr = headfits(filen)
dtline=where(strpos(hdr,keyword) EQ 0)
; if there are two instances of the keyword, choose the first
if(n_elements(dtline) GT 1) then dtline = dtline[0]

stspec = hdr[dtline]
; get the date string and make it into UT
tspec = strsplit(stspec,"'",/extract)

JDspec = date_conv(tspec[1],'J')

return, JDspec

end
