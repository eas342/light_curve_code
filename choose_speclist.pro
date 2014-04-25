pro choose_speclist
;; Chooses which spectrum to look at and saves to current_speclist.txt

fchoice = choose_file(searchDir='file_lists',filetype='.txt')
spawn,'cp '+fchoice+' file_lists/current_speclist.txt'

;; If it's KIC 1255 data, search for an associated date to
;; recall the MORIS photometry
startpos = strpos(fchoice,'kic1255') + 8
if startpos NE -1 then begin
   useDate = strmid(fchoice,startpos,9)
endif else useDate = 'NA'
SpecListName = fchoice
save,useDate,SpecListName,filename='data/used_date.sav'

end
