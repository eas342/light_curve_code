pro timeser,showimg=showimg,inl=inl,boxfn=boxfn,doflat=doflat
;; generate a time series for a given region

;showimg shows the spectra files
;inl - input spectral list could be one of the following:
;'cut1' skips the first image in each set because it always appears lower
;;'boxfn' allows you to choose which region file to choose for box regions
;; dowav keyword will save wavelengths for the boxes

common commonV,nbox

startc1 = systime(/seconds) ;; start the clock
;; check if the inl parameter is undefined
if n_elements(inl) EQ 0 then inl = 'file_lists/low_lamp_speclist.txt'
if n_elements(boxfn) EQ 0 then boxfn = 'regions/start_boxes_ord4-6.reg'

print,"Analyzing spectra listed in ",inl,"..."      
readcol,inl,filen,format='(A)',/silent,comment='#' ;; input file for the list of spectra to look at
;; files with a comment symbol are ignored


;; check if the boxchoice parameter is undefined
;if n_elements(boxch) EQ 0 then boxch = 'orig'
;case boxch of
;   'HD8991A': boxfn = 'K-bandHD8991A.reg'
;   else: begin
;      boxfn = 'K-bandHD8991A.reg'
;      boxfnB = 'K-bandHD8991B.reg'
;   end
;endcase
readcol,boxfn,ctype,rtype,xcent,ycent,xsize,ysize,zeros,$
        format='(A,A,F,F,L,L,L)',/silent ;; input regions for flux extractions

nbox = n_elements(xcent)

nfile = n_elements(filen)
psymarr = [1,2,4,5,6,7] ;; choice of 6 symbol times
t = dblarr(nfile) ;; for the time in JD
flarr = dblarr(nfile,nbox) ;; array to keep all the fluxes for given boxes
totfl = dblarr(nfile,nbox) ;; array to keep track of the total flux
flerr = dblarr(nfile,nbox) ;; array for the photon errors in flux
ipx = lonarr(1024,1024l) ;; array to keep track of the pixels used by the boxes
expt = dblarr(nfile) ;; for the exposure times
wavl = dblarr(nbox) ;; for the wavelengths



if keyword_set(doflat) then begin
   flatnm = '../../corot1/IRTF_UT2011Dec23_corot1_bigdog_guidedog/raw/bigdog/cal/masterflat1.fits'
   ;; read the flat field (which is Y-inverted by SpeXTool)
   flatInv = mrdfits(flatnm,/dscale)
   ;; flip the flat over the y axis (b/c SpeXTool flipped it)
   flat = flatInv[*,reverse(lindgen(1024))]

   ;; normalize
   flat = flat / mean(flat)
   nonzeropx = where(flat NE 0) ;; non zero pixels
   zeropx = where(flat EQ 0)    ;; zero pixels

   assert,(n_elements(nonzeropx)+n_elements(zeropx)),'=',1024l*1024l,'All pixels not accounted for'

   if keyword_set(showimg) then begin

      expand,flat,512,512,ims ;; ironically, to shrink
      
      st = sort(ims)
      maxv = ims[st[round(0.99*n_elements(st))]]
      minv = ims[st[round(0.05*n_elements(st))]]
      window,0,xsize=512,ysize=512
      
      tv,hist_equal(ims*255l/(maxv-minv)),0,0 ; switch to bytes for display purposes
      imfilen = 'plots/img_flat.png'
      write_png,imfilen,tvrd(true=1)
   endif

endif


for i=0l,nfile - 1l do begin
   print,'Extracting Boxes for image ',i+1,' of ',nfile
   img = mrdfits(filen[i],/silent,0,hdr,/dscale)
   npx = n_elements(img)
   ;; flat field
   if keyword_set(doflat) then begin
      img[nonzeropx] = img[nonzeropx] / flat[nonzeropx]
      img[zeropx] = 0
   endif

   ;; get the MJD
   t[i] = fxpar(hdr,'MJD_OBS')
;   t[i] = date_conv(MJD,'FITS')
   hdr = headfits(filen[i])
   expt[i] = fxpar(hdr,'ITIME')


   if i EQ 0l then begin ;; only figure out coordinates the first time
      coord = array_indices(img,lindgen(npx))
      xcoor = coord[0,*]
      ycoor = coord[1,*]

      ;; Define extraction boxes
      for j=0l,nbox-1l do begin
              xcentC = xcent
              ycentC = ycent
              xsizeC = xsize
              ysizeC = ysize
         XY = boxcoor(xcentC[j],ycentC[j],sizeX=xsizeC[j],sizeY=ysizeC[j])

         bpx = where(xcoor GE XY[0,0] and xcoor LT XY[0,1] and $
                     ycoor GE XY[1,0] and ycoor LT XY[1,1])

         ipx[bpx] = j + 1
         ;; save the box pixels to a structure
         if j EQ 0l then begin
            sbpx = create_struct('bpx0',bpx)
         endif else sbpx = create_struct(sbpx,'bpx'+strtrim(j,1),bpx)
 
         ;; Find wavelengths for the boxes
;         wavl[j] = boxcoor(xcent[j],ycent[j],sizeX=xsize[j],sizeY=ysize[j],/wavelength)
;         print,'Box ',j,' wavelength = ',wavl[j]
      endfor
      forprint,textout='data/bx_wavelengths.txt',lindgen(nbox),wavl,format='(A,F11.4)',comment="# Box#   Wavelength (microns)",/silent
   endif
   ;; Extract the fluxes and save them to file
   for j=0l,nbox-1l do begin
      flarr(i,j) = total(img[sbpx.(j)]) ;; image i and box number j
      ;; keep the total value for use in error estimation
      totfl(i,j) = flarr(i,j)
      ;; divide by the number of pixels
      flarr(i,j) = flarr(i,j) / float(n_elements(img[sbpx.(j)]))
   endfor


   if keyword_set(showimg) and (i EQ 0l ) then begin
      red = 255l
      green = 255l*255l
      lblue = 256l*256l*255l+255l*256l
      white = 256l*256l*255l+256l*255l+255l
      magenta = 256l*256l*255l+255l
      yellow = 255l+255l*256l

      colorOpt = [white,red,green,lblue,magenta,yellow]
      colorarr = colorOpt[lindgen(nbox) mod 6]
      lineOpt = [0,2,3,4,5]
      linarr = lineOpt[lindgen(nbox) mod 4]


      ;; Blank out the boxes
;      boxp = where(ipx) GT 0l
;      img[boxp] = 0 ;; zero out the boxes


      ;; Prep the image for viewing
      ;; shrink the image (ironically with the expand command!) for viewing purposes
      expand,img,512,512,ims
      
      ;; also make an image of the boxes
      expand,ipx,512,512,bpxs

      st = sort(ims)
      maxv = ims[st[round(0.99*n_elements(st))]]
      minv = ims[st[round(0.05*n_elements(st))]]
      window,0,xsize=512,ysize=512
      
      tv,hist_equal(ims*255l/(maxv-minv)),0,0 ; switch to bytes for display purposes
;      tv,(ims*255l/(maxv-minv)),0,0 ; switch to bytes for display purposes
;      tv,hist_equal(bpxs),0,0 ;show boxes
      for j=0l,nbox-1l do begin

         ;; Draw the boxes
         XY = boxcoor(xcentC[j],ycentC[j],sizeX=xsizeC[j],sizeY=ysizeC[j],/fordraw)
;         XY = boxcoor(xcent[j],ycent[j],sizeX=xsize[j],sizeY=ysize[j],/fordraw)
         plots,XY[0,*],XY[1,*],color=0,/device,thick=3.7
         plots,XY[0,*],XY[1,*],color=colorArr[j],/device,$
               linestyle=linArr[j],thick=2.5

         ;; Label the boxes
;         oplot,[xcent[j]/2d,ycent[j]/2d],psym=
         xyouts,xcent[j]/2d,ycent[j]/2d,strtrim(j,1),/device,charsize=1.6,$
                color=white,charthick=3
         xyouts,xcent[j]/2d,ycent[j]/2d,strtrim(j,1),/device,charsize=1.6,$
                color=red

      endfor
      imfilen = 'plots/box_choices.png'
      write_png,imfilen,tvrd(true=1)

   endif

endfor

;; calculate the expected photon errors
;Gain = float(fxpar(hdr,'GAIN'))
Gain = 13.0 ;; looked up on Rayner et al. 2003
errfrac = 1E / sqrt(totfl * Gain) ;sqrt(N) / N (=fractional error)

;; multiply the fractional errors by the counts per pixel
flerr = errfrac * flarr

;; save the fluxes
openw,1,'data/fl_time_ser.txt'
printf,1,'# Time (JDB)      Exposure Time     average DN/s'
for i=0l,nfile-1l do begin
   printf,1,t[i],format='($,F23.12)' ;; save the start time
   printf,1,expt[i],format='($,F15.6)' ;; save the exposure time
   for j=0l,nbox-1l do begin
      ; make sure to add a newline where needed
      if j EQ nbox-1l then fmk='(F22.2)' else fmk='($,F22.2)'
      printf,1,flarr(i,j),format=fmk
   endfor
endfor
close,1

;; save the expt,time and flux arrays, and erros
save,filename='data/fluxsave.sav',flarr,t,expt,nbox,nfile,flerr,wavl,t

;; save the photon errors for the fluxes
;; Find the photon errors

end

