pro prep_imcentroid_for_geomap
  ;; get the number of files
  readcol,'../IRTF_UT2012Jan04/file_lists/guider100list.txt',origname,format='(A)'
  nfile = n_elements(origname)

  ;; get the number of sources
  readcol,'../IRTF_UT2012Jan04/file_lists/astrom_starlist.txt',firstcol,format='(A)',skipline=2
  nstars = n_elements(firstcol)

  ;; get the reference data
  imcentroidNM = '../IRTF_UT2012Jan04/phot_data/centroids100.txt'
  centroidDIR = '../IRTF_UT2012Jan04/phot_data/file_by_file_cen'

  readcol,imcentroidNM,filenfull,xref,xEref,yref,yEref,starnumref,$
          format='(A,F,A,F,A,I)',skipline=1l,$
          numline=nstars

  ;; Read in the centroids one file at a time
  for i=0l,nfile-1l do begin
     ;; get the centroids
     readcol,imcentroidNM,filenfull,x,xE,y,yE,starnum,$
             format='(A,F,A,F,A,I)',skipline=1l+(nstars+1l)*i,$
             numline=nstars,/silent
     filenarr = strsplit(filenfull[0],'/',/extract)
     nterms = n_elements(filenarr)

     ;; make file
     openw,1,centroidDIR+'/centroids_'+filenarr[nterms-1]+'_.dat'
     for j=0l,nstars-1l do begin
        printf,1,xref[j],' ',yref[j],' ',x[j],' ',y[j],$
               format='(F,A,F,A,F,A,F)'
     endfor
     close,1
  endfor
end
