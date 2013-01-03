pro ev_flatfield
;; Makes a flat field from a set of images

  readcol,'file_lists/flat_list.txt',flatlistnm,format='(A)'
  nfiles = n_elements(flatlistnm)

  ;; get the image size from the first file
  firstHeadr = headfits(flatlistnm[0])

  Xlength = sxpar(firstHeadr,'NAXIS1')
  Ylength = sxpar(firstHeadr,'NAXIS2')
  
  ;; Create an array for the flat field
  flatarray = fltarr(Xlength,Ylength)
  
  ;; Create an array for all the data
  dataarray = lonarr(Xlength,Ylength,nfiles)

  ;; read in each file
  for i=0l,nfiles-1l do begin
     img = readfits(flatlistnm[i],/silent)
     dataarray[*,*,i] = img
  endfor
  
  ;; median combine
  flatarray = median(dataarray,dimension=3)

  ;; normalize
  flatarray = flatarray / float(median(flatarray))

  ;; display
  window,1,xsize=Xlength,ysize=Ylength
  ;; find smallest 0.5% of values
  ;; find largest 99.5% of values
  sortpoints = sort(flatarray)
  sortflat = flatarray[sortpoints]
  mindisplay = sortflat[Xlength*Ylength*5l/1000l]
;  maxdisplay = sortflat[Xlength*Ylength*995l/1000l]
  maxdisplay = sortflat[Xlength*Ylength*600l/1000l]
  tvscl,(flatarray-mindisplay)/maxdisplay*255b
  writefits,'cals/median_flat.fits',flatarray
  stop
end
