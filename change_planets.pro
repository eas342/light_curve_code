pro change_planets

  readcol,'data/planet_list.txt',planets,format='(A)'
  nplanets = n_elements(planets)
  for i=0l,nplanets-1l do begin
     print,string(i,format='(I5)'),' ',planets[i]
  endfor
  read,'Planet choice: ',planetchoice

  pname = planets[planetchoice]

  spawn,'cp transit_info/transit_epoch_'+pname+'.txt transit_info/transit_epoch.txt'
  spawn,'cp transit_info/planet_info_'+pname+'.txt transit_info/planet_info.txt'
  spawn,'cp data/object_coordinates_'+pname+'.txt data/object_coordinates.txt'
  spawn,'cp param_input/specphot_range_'+pname+'.txt param_input/specphot_range.txt'
  spawn,'cp param_input/time_series_sep_'+pname+'.txt param_input/time_series_sep.txt'

end
