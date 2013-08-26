pro get_profile_widths,showplot=showplot
;; Fits Gaussians to to the spatial profeils in the spectrum
;; showplot -- show a plot, otherwise it just records the numbers

  readcol,'file_lists/current_speclist.txt',fileL,format='(A)',/silent

  nfile = n_elements(fileL)

  Widths = fltarr(nfile,2)

  for j=0l,nfile-1l do begin

     ;; Find the original file name
     endS = strpos(fileL[j],'.ms.d.fits')
     origNm = strmid(fileL[j],0,endS)+'_straight.fits'
     print,origNm
     a = mrdfits(origNm,0,header,/silent)
     b = total(a,1)

     indices=findgen(n_elements(b))
     profstruct = create_struct('data',[[indices],[b]])
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

        if keyword_set(showplot) then begin
           plot,indices,b
           oplot,subarrayX,result,color=mycol('yellow'),thick=2
;     oplot,subarrayX,A[0] * exp(-((subarrayX - A[1])/A[2])^2) + A[3]
;                                               + A[4] * SubarrayX,thick=2,color=mycol('yellow')
        endif
     endfor
  endfor
  save,Widths,filename='data/prof_widths.sav'
end

