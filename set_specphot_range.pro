pro set_specphot_range,lowF,highF
  openw,1,'param_input/specphot_range.txt'
  printf,1,'# Range Start   Range End'
  printf,1,lowF,highF,format='(2F10.3)'
  close,1
end
