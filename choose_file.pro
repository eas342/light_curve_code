function choose_file,searchDir=searchDir,filetype=filetype
;; Returns a user-chosen file from a directory. The default directory
;; is the current one IDL was run from

  ;; get the current directory
  cd,c=currentd
  ;; Search the given directory
  if n_elements(searchDir) EQ 0 then searchDir = '' else begin
     searchDir = searchDir+'/'
  endelse
  if n_elements(filetype) EQ 0 then filetype = ''

  fileopt = file_search(currentd+'/'+searchDir+'*'+filetype)
  print,'Available Files. Choose one of the numbers:'
  for i=0l,n_elements(fileopt)-1l do begin
     trimst = strsplit(fileopt[i],'/',/extract)
     print,string(i,Format='(I5)'),' ',trimst(n_elements(trimst)-1l)
  endfor
  read,'File Choice: ',filechoice

  return, fileopt[filechoice]
end
