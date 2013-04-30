pro straighten_spec,inlist,outlist

;; Read in the straightening data
readcol,'data/shift_data/shift_vals_from_arc.txt',rowN,shiftMod,format='(F,F)',skipline=1

;; Read the list of spectra to straighten & their output file names
readcol,inlist,infiles,format='(A)'
readcol,outlist,outfiles,format='(A)'
nfiles = n_elements(infiles)
assert,nfiles,'=',n_elements(outfiles),'In/Out file lists are mismatched.'

for j=0l,nfiles-1l do begin
;   img =
;   mrdfits('../IRTF_UT2012Jan04/proc/bigdog/bigdog0001.a.fits',0,origHeader)
   img = mrdfits(infiles[j],0,origHeader)

   imgSize = size(img)
   NX = imgSize[1]
   NY = imgSize[2]
   assert,n_elements(rowN),'=',NY,"Shift row and image row numbers don't match"
   recimg = fltarr(NX,NY)
   
   for i=0l,NY-1l do begin
      recimg[*,i] = shift_interp(img[*,i],shiftMod[i])
   endfor
   
   ;; Check if file exists so you don't overwrite anything!
   checkfile = findfile(outfiles[j],count=count)
   if count GE 1 then begin
      print,"Out File Exists!"
      stop
   endif else writefits,outfiles[j],recimg,origHeader
endfor

end
