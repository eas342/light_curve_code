function ev_delim_read,filen,delimiter=delimiter,only1header=only1header
; reads any data file with data separated by a delimiter
; data is returned as a structure where the tag names are the first
; line of the data file
; for a CSV file, the delimeter should be set to ','
; filen is the file name one wishes to reach
; delimiter is the data delimeter - default is pipe '|'
; only1header -- only one header line (no long,double,string descriptions)

if n_elements(delimiter) EQ 0 then delimiter = '|'
nrows = file_lines(filen)
;; get the tag names
openr,1,filen
tagst = ''
readf,1,tagst
tags = strtrim(strsplit(tagst,delimiter,/extract),1)
ntags = n_elements(tags)

replaceTypes = [' ','(',')','/','\','.','-']
nTypes = n_elements(replaceTypes)
for i=0l,ntags-1l do begin
   ;; Trim out the parentheses and spaces, etc
   for j=0l,nTypes-1l do begin
      tags[i] = strjoin(strsplit(tags[i],replaceTypes[j],/extract),'_')
   endfor
;   print,tags[i]
endfor

if not keyword_set(only1header) then begin
typest = ''
readf,1,typest
typest = strtrim(strsplit(typest,delimiter,/extract),1)
;; Change all funky variable types into strings
funnyp = where(typest EQ 'numeric' or typest EQ 'datetime')
if funnyp NE [-1] then typest[funnyp] = 'string'
funnyp = where(typest EQ 'ra' or typest EQ 'dec')
if funnyp NE [-1] then typest[funnyp] = 'double'
funnyp = where(typest EQ 'integer')
if funnyp NE [-1] then typest[funnyp] = 'long'
endif

;; get the data
filed = ''
file2d = strarr(nrows-2l,ntags)
for i=0l,nrows-3l do begin
   readf,1,filed
   file2d[i,*] = strtrim(strsplit(filed,delimiter,/extract,/preserve_null),1)
endfor

close,1
free_lun,1

;; make structure
st = create_struct(tags[0],file2d[*,0])
for j=1l,ntags-1l do begin
   ;; convert to the correct file type
   junk = execute('ColumnArr = '+typest[j]+'(file2d[*,j])')
   st = create_struct(st,tags[j],columnArr)
endfor

return,st
end
