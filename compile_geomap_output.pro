pro compile_geomap_output
  ;; get the number of files
;  readcol,'../IRTF_UT2012Jan04/file_lists/proc100list.txt',origname,format='(A)'
  readcol,'../IRTF_UT2012Jan04/file_lists/fullguider_list.txt',origname,format='(A)'
  nfile = n_elements(origname)

  ;; Set up data arrays
  xshift = dblarr(nfile)
  xrms = dblarr(nfile)
  yshift = dblarr(nfile)
  yrms = dblarr(nfile)
  xrot = dblarr(nfile) ;; x Rotation in degrees
  yrot = dblarr(nfile)
  filenl = strarr(nfile) ;; list of file names

  ;; get the number of sources
;  readcol,'../IRTF_UT2012Jan04/file_lists/astrom_starlist.txt',firstcol,format='(A)',skipline=2
;  nstars = n_elements(firstcol)

  datalength = 28

  ;; get the reference data
;  geomapNM = '../IRTF_UT2012Jan04/phot_data/geomap_out.dat'
  geomapNM = '../IRTF_UT2012Jan04/phot_data/geomap_outfull.dat'
  outputNM = '../IRTF_UT2012Jan04/phot_data/geomap_listfull.txt'

  ;; Read in the centroids one file at a time
  for i=0l,nfile-1l do begin
     ;; get the output for that file
     readcol,geomapNM,description,value,$
             format='(A,A)',skipline=1l+(datalength+1l)*i,$
             numline=datalength,/silent
     filenarr = strsplit(value[0],'/',/extract)
     nterms = n_elements(filenarr)
     filenl[i] = filenarr[nterms-1]
     xshift[i] = double(value[7])
     yshift[i] = double(value[8])
     xrot[i] = double(value[11])
     yrot[i] = double(value[12])
     xrms[i] = double(value[13])
     yrms[i] = double(value[14])
  endfor

  ;; for the rotations, make angles more than 180 degrees negative
  highAng = where(xrot GT 180D)
  if highAng NE [-1] then xrot[highAng] = xrot[highAng] - 360D
  highAng = where(yrot GT 180D)
  if highAng NE [-1] then yrot[highAng] = yrot[highAng] - 360D
  

  ;; Save the data
  forprint,textout=outputNM,$
           comment="#Filename,  X Shift(px), X RMS,  Y Shift, Y RMS,  X Rot(deg),  Y Rot,",$
           filenl,xshift,xrms,yshift,yrms,xrot,yrot,$
           width=200
  
  stop
end
