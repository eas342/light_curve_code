pro choose_speclist,fchoice=fchoice
;; Chooses which spectrum to look at and saves to current_speclist.txt

if n_elements(fchoice) EQ 0 then begin
   fchoice = choose_file(searchDir='file_lists',filetype='.txt')
endif
spawn,'cp '+fchoice+' file_lists/current_speclist.txt'

;; If it's KIC 1255 data, search for an associated date to
;; recall the MORIS photometry

case 1 of
   strpos(fchoice,'kic1255') NE -1: begin
      startpos = strpos(fchoice,'kic1255') + 8
      useDate = strmid(fchoice,startpos,9)
      dateLength = strlen(useDate)
      showDate = strmid(useDate,0,4) + ' ' + strupcase(strmid(useDate,4,1)) +$
                 strmid(useDate,5,dateLength - 5)
   end
   strpos(fchoice,'corot2') NE -1: begin
      startpos = strpos(fchoice,'corot2') + 7
      usedate = strmid(fchoice,startpos,10)
      showdate = strjoin(strsplit(usedate,'_',/extract),'-')

   end
   else: begin
      useDate = 'NA'
      showDate = 'NA'
   end
endcase

SpecListName = fchoice

splitName = strsplit(speclistname,'/',/extract)
nsplitNames = n_elements(splitName)
splitFileName = strsplit(splitName[nsplitNames-1],'.',/extract)
specfileListNamePrefix = splitFileName[0]

save,useDate,SpecListName,specfileListNamePrefix,showDate,$
     filename='data/used_date.sav'

end
