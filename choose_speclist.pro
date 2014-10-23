pro choose_speclist,fchoice=fchoice
;; Chooses which spectrum to look at and saves to current_speclist.txt

if n_elements(fchoice) EQ 0 then begin
   fchoice = choose_file(searchDir='file_lists',filetype='.txt')
endif
spawn,'cp '+fchoice+' file_lists/current_speclist.txt'

;; If it's KIC 1255 data, search for an associated date to
;; recall the MORIS photometry
startpos = strpos(fchoice,'kic1255') + 8
if startpos NE -1 then begin
   useDate = strmid(fchoice,startpos,9)
endif else useDate = 'NA'
SpecListName = fchoice

splitName = strsplit(speclistname,'/',/extract)
nsplitNames = n_elements(splitName)
splitFileName = strsplit(splitName[nsplitNames-1],'.',/extract)
specfileListNamePrefix = splitFileName[0]

save,useDate,SpecListName,specfileListNamePrefix,filename='data/used_date.sav'

end
