pro choose_speclist
;; Chooses which spectrum to look at and saves to current_speclist.txt

fchoice = choose_file(searchDir='file_lists',filetype='.txt')
spawn,'cp '+fchoice+' file_lists/current_speclist.txt'

end
