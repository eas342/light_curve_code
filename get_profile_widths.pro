pro get_profile_widths,showplot=showplot,jan04corot1=jan04corot1,$
                       dec23corot1=dec23corot1,dec29corot1=dec29corot1
;; Fits Gaussians to to the spatial profeils in the spectrum
;; showplot -- show a plot, otherwise it just records the numbers
;; jan04/dec23/dec29corot1 -- works for Corot-1 where I have different file name conventions


  readcol,'file_lists/current_speclist.txt',fileL,format='(A)',/silent,$
          stringskip='#'

  nfile = n_elements(fileL)

  Widths = fltarr(nfile,2);; planet host star =0, reference star=1
  starLocations = fltarr(nfile,2)

  case 1 of 
     keyword_set(jan04corot1): begin
        searchKey = '.ms.fits'
        imageName = '.a.fits'
     end
     keyword_set(dec23corot1): begin
        searchKey = '.d.fits'
        imageName = '.a.fits'
     end
     keyword_set(dec29corot1): begin
        searchKey = '.a.ms.d.fits'
        imageName = '.a.fits'
     end
     else: begin
        searchKey = '.ms.d.fits'
        imageName = '_straight.fits'
     end
  endcase

  for j=0l,nfile-1l do begin

     ;; Find the original file name
     endS = strpos(fileL[j],searchKey)
     origNm = strmid(fileL[j],0,endS)+imageName

     print,origNm
     a = mrdfits(origNm,0,header,/silent)
     b = total(a,1)

     indices=findgen(n_elements(b))
     ;; smooth the profile first
     smoothprof = smooth(b,10)
     profstruct = create_struct('data',[[indices],[smoothprof]])
     mc_findpeaks,profstruct,2,1,posit,apsign,/auto

;     print,posit
     HW = 40 ;; Half width
     for i=0l,1 do begin
        ;; Go through both sources
        ;; make a sub-array surrounding a peak
        startx = posit[i] - HW
        if startx LT 0l then startx = 0l
        endx = posit[i] + HW
        if endx GE n_elements(b) -1l then endx = n_elements(b) - 1l
        subarrayX = indices[startx:endx]
        subarrayY = b[startx:endx]
        result = gaussfit(subarrayX,subarrayY,A,nterms=5)
        Widths[j,i] = A[2]
        starLocations[j,i] = A[1]

        if keyword_set(showplot) then begin
           plot,indices,b
           oplot,subarrayX,result,color=mycol('yellow'),thick=2
;     oplot,subarrayX,A[0] * exp(-((subarrayX - A[1])/A[2])^2) + A[3]
;                                               + A[4] *
;                                               SubarrayX,thick=2,color=mycol('yellow')
           ;; stop
           if j EQ 50 then stop
        endif
     endfor
  endfor
  save,Widths,starLocations,filename='data/prof_widths.sav'
end

