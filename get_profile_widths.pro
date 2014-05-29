pro get_profile_widths,showplot=showplot,jan04corot1=jan04corot1,$
                       dec23corot1=dec23corot1,dec29corot1=dec29corot1,$
                       useSaved=useSaved,autofind=autofind
;; Fits Gaussians to to the spatial profeils in the spectrum
;; showplot -- show a plot, otherwise it just records the numbers
;; jan04/dec23/dec29corot1 -- works for Corot-1 where I have different file name conventions
;; usedSaved -- used saved data for a given speclist
;; autofind - find the peaks (instead of using the aperture centers as
;;            starts)

  readcol,'file_lists/current_speclist.txt',fileL,format='(A)',/silent,$
          stringskip='#'
  ;; also get the speclist name
  restore,'data/used_date.sav'

  nfile = n_elements(fileL)

  if not keyword_set(useSaved) then begin
     Widths = fltarr(nfile,2) ;; planet host star =0, reference star=1
     widthsE = widths
     starLocations = fltarr(nfile,2)
     starLocationsE = starlocations
     voigts = fltarr(nfile,2)
     voigtsE = voigts

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
        endS = strpos(fileL[j],searchKey) ;; where the file name ends
        startS = strpos(fileL[j],'/',/reverse_search) ;; where the file name starts
        origNm = strmid(fileL[j],0,endS)+imageName ;; original 2d image name + directory

        ;; image name stuff without extension (usually _straight)
        postpos = strpos(imageName,'.')
        ;; Name of the aperture file
        apFileNm = strmid(fileL[j],0,startS+1)+'database/ap'+$
                   strmid(fileL[j],startS+1,(endS-startS-1))+strmid(imageName,0,postpos)
        
        
        if j mod 20 EQ 0 then print,origNm
        a = mrdfits(origNm,0,header,/silent)
        b = total(a,1)
        
        indices=findgen(n_elements(b))
        ;; smooth the profile first
        smoothprof = smooth(b,10)
        profstruct = create_struct('data',[[indices],[smoothprof]])

        if keyword_set(autofind) then begin
           mc_findpeaks,profstruct,2,1,posit,apsign,/auto
        endif else begin
           ;; Read from aperture file
           line = ''
           i=0
           posit = fltarr(2)
           openr,1,apFilenm
           while not eof(1) do begin
              readf,1,apFilenm,line
              if strpos(line,'center') NE -1 then begin
                 centerArr = strsplit(line,' ',/extract)
                 posit[i] = float(centerArr[1])
                 i = i+1
              endif
           endwhile
           free_lun,1
        endelse
        HW = 40 ;; Half width

        if abs(posit[1] - posit[0]) LT HW then begin
           ;; If the two stars are close, we'll need to fit
           ;; them simultaneously instead of one at a time
           if max(posit) GT max(indices) or min(posit) LT 0 then begin
              print,'Peak finding Failed in Image',origNm
              starLocations[j,i] = !values.f_nan
              Widths[j,i] = !values.f_nan
           endif
           
           HW = 80
           startx = min(posit) - HW
           if startX LT 0l then startx=0l
           endx = max(posit) + HW
           if endx GE n_elements(b) -1l then endx = n_elements(b) - 1l

           ;; Estimate outside points for error
           outp = where(indices LT startx or indices GT endx)
           yerr = robust_sigma(b[outp]) + fltarr(endx-startx+1)
           subarrayX = indices[startx:endx]
           subarrayY = b[startx:endx]
;           expr = '(voigt_approx(P[0],X - P[1],P[2]) * P[6]+'+$
;                  ' P[7] * voigt_approx(P[3],X - P[4],P[5]))'+$
;                  '+ eval_legendre(X,P[8:10])'
           expr = '(voigt_approx(P[0],X - P[1],P[2]) * P[6]+'+$
                  ' P[7] * voigt_approx(P[0],X - P[4],P[2]))'+$
                  '+ eval_legendre(X,P[8:10])'
;           start = [0.1E,posit[0],3E,0.1E,posit[1],3E,max(b) *
;           3E,max(b) * 3E,max(b) * 3E,0E,0E]
           start = [0.1E,posit[0],3E,0.1E,posit[1],3E,max(b) * 3E,max(b) * 3E,max(b) * 3E,0E,0E]
           nparams = n_elements(start)
           pi = replicate({fixed:0, limited:[0,0], limits:[0.0E,0.0E]},nparams)
;           fixedP = [0]
;           pi[fixedP].fixed = 1
           pi[0:5].limited = [1,1] ;; limit the profile parameters
           indV = [0E,3E] ;; voigt parameter Indices
           pi[indV].limits = [0E,5E]
           indPos = [1,4] ;; Position indices
           pi[indPos].limits = [-startx,endx]
           indSig = [2,5] ;; Sigma indices
           pi[indSig].limits = [-HW/2E,HW/2E]
           
           result = mpfitexpr(expr,subarrayX,subarrayY,yerr,start,parinfo=pi,perr=perr,/quiet)
           yfit = expression_eval(expr,subarrayX,result)
           resid = yfit - b
           
           if keyword_set(showplot) then begin
              plot,indices,b,xtitle='Y position (px)',$
                   ytitle='Flux'
              oplot,subarrayX,yfit,color=mycol('yellow'),thick=2
           endif
           Widths[j,0] = result[2]
           WidthsE[j,0] = perr[2]
           Widths[j,1] = result[2];result[5]
           WidthsE[j,1] = perr[2];perr[5]
           starlocations[j,0] = result[1]
           starlocationsE[j,0] = perr[1]
           starlocations[j,1] = result[4]
           starlocationsE[j,1] = perr[4]
           voigts[j,0] = result[0]
           voigtsE[j,0] = perr[0]
           voigts[j,1] = result[0];result[3]
           voigtsE[j,1] = perr[0];perr[3]

           acheck = widths[0:j,0]
;           if j GT 20 and $
;              acheck[j] - median(acheck) GT 10E *
;                          robust_sigma(acheck)then stop
;           if j GT 30 then stop
        endif else begin
           for i=0l,1 do begin
              ;; Go through both sources
              ;; make a sub-array surrounding a peak
              if posit[i] GT max(indices) or posit[i] LT 0 then begin
                 print,'Peak finding Failed in Image',origNm
                 starLocations[j,i] = !values.f_nan
                 Widths[j,i] = !values.f_nan
              endif else begin
                 startx = posit[i] - HW
                 if startx LT 0l then startx = 0l
                 endx = posit[i] + HW
                 if endx GE n_elements(b) -1l then endx = n_elements(b) - 1l
                 subarrayX = indices[startx:endx]
                 subarrayY = b[startx:endx]
                 result = gaussfit(subarrayX,subarrayY,A,nterms=5)
                 Widths[j,i] = A[2]
                 starLocations[j,i] = A[1]
              endelse
              if keyword_set(showplot) then begin
                 if i EQ 0 then begin
                    zsubarrayX = subarrayX
                    zresult = result
                 endif else begin
                    plot,indices,b
                    oplot,subarrayX,result,color=mycol('yellow'),thick=2
;     oplot,subarrayX,A[0] * exp(-((subarrayX - A[1])/A[2])^2) + A[3]
;                                               + A[4] *
;                                               SubarrayX,thick=2,color=mycol('yellow')
                    ;; stop
                    oplot,zsubarrayX,zresult,color=mycol('lblue'),thick=2
                    if j GE 50 then stop
                 endelse
              endif
           endfor
        endelse
        
     endfor
     save,Widths,starLocations,voigts,$
          widthsE,starlocationsE,voigtsE,$
          filename='data/state_parameters/widths/widths_'+$
          specfileListNamePrefix+'.sav'
  endif

  restore,'data/state_parameters/widths/widths_'+$
          specfileListNamePrefix+'.sav'
  save,Widths,starLocations,voigts,$
       widthsE,starlocationsE,voigtsE,$
       filename='data/prof_widths.sav'

end


