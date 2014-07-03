pro get_profile_widths,showplot=showplot,jan04corot1=jan04corot1,$
                       dec23corot1=dec23corot1,dec29corot1=dec29corot1,$
                       useSaved=useSaved,autofind=autofind,$
                       specific=specific,psplot=psplot
;; Fits Gaussians to to the spatial profeils in the spectrum
;; showplot -- show a plot, otherwise it just records the numbers
;; jan04/dec23/dec29corot1 -- works for Corot-1 where I have different file name conventions
;; usedSaved -- used saved data for a given speclist
;; autofind - find the peaks (instead of using the aperture centers as
;;            starts)
;; specific - choose a specific file to fit profiles of (usually use
;;            with showplot
;; psplot - save a postscript plot

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/spatial_profiles/spatial_prof_'+$
                 string(specific,format='(I03)')
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
     device,xsize=15, ysize=10,decomposed=1,/color
     !p.thick=3.5
     !x.thick=2
     !y.thick=2
     thickline = 5
     thinline = 3.5
  endif else begin
     thickline=2
     thinline = 1
  endelse

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
     starAmps = fltarr(nfile,2) ;; amplitudes of the stars
     starAmpsE = fltarr(nfile,2) 

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
 
     if n_elements(specific) EQ 0 then begin
        startFile=0l
        endFile=nfile-1l
     endif else begin
        startFile = specific
        endfile=startfile
     endelse 
     for j=startFile,endFile do begin
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
        b = total(a[60:159,*],1)

        b = b/max(b)

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
        HW = 40 ;; Half width for looking at individual stars

        if abs(posit[1] - posit[0]) LT HW then begin
           ;; If the two stars are close, we'll need to fit
           ;; them simultaneously instead of one at a time
           if max(posit) GT max(indices) or min(posit) LT 0 then begin
              print,'Peak finding Failed in Image',origNm
              starLocations[j,i] = !values.f_nan
              Widths[j,i] = !values.f_nan
           endif
           
           HWBoth = 80 ;; half width of baseline for two stars
           startx = min(posit) - HWboth
           if startX LT 0l then startx=0l
           endx = max(posit) + HWboth
           if endx GE n_elements(b) -1l then endx = n_elements(b) - 1l

           ;; Estimate outside points for error
           outp = where(indices LT startx or indices GT endx)
           yerr = robust_sigma(b[outp]) + fltarr(endx-startx+1)
           subarrayX = indices[startx:endx]
           subarrayY = b[startx:endx]
           expr = '(voigt_approx(P[0],X - P[1],P[2]) * P[6]+'+$
                  ' P[7] * voigt_approx(P[3],X - P[4],P[5]))'+$
                  '+ eval_legendre(X,P[8:10])'
           start = [0.1E,posit[0],3E,0.1E,posit[1],3E,max(b) * 3E,max(b) * 3E,max(b) * 3E,0E,0E]
           nparams = n_elements(start)
           pi = replicate({fixed:0, limited:[0,0], limits:[0.0E,0.0E]},nparams)
           indV = [0E,3E] ;; voigt parameter Indices
           pi[indV].limited = [1,0]
           pi[indV].limits = [0E,0E]
           indPos = [1,4] ;; Position indices
           pi[indPos].limited = [1,1]
           pi[indPos].limits = [-startx,endx]
           indSig = [2,5] ;; Sigma indices
           pi[indSig].limited = [1,1]
           pi[indSig].limits = [-HWboth/2E,HWboth/2E]
           
           result = mpfitexpr(expr,subarrayX,subarrayY,yerr,start,parinfo=pi,perr=perr,/quiet)
           yfit = expression_eval(expr,subarrayX,result)
           resid = yfit - b[startx:endx]
           
           if keyword_set(showplot) then begin
              !p.multi = [0,1,2]
              plot,subarrayX,subarrayY,xtitle='Y position (px)',$
                   ytitle='Flux'
              oplot,subarrayX,yfit,color=mycol('blue'),thick=thickline

              al_legend,['Measured','Model'],linestyle=[0,0],$
                        color=[!p.color,mycol('blue')],$
                        thick=[thinline,thickline]
              al_legend,[cgGreek('sigma')+' = '+$
                         string(result[2],format='(F6.2)')+' +/- '+$
                         string(perr[2],format='(F6.2)'),$
                         'a = '+$
                         string(result[0],format='(F6.2)')+' +/- '+$
                         string(perr[0],format=('(F6.2)'))],$
                        /right
;                         'a = '+string(result[
              plot,subarrayX,subarrayY,xtitle='Y position (px)',$
                   ytitle='Residual',yrange=threshold(resid),/nodata
              oplot,subarrayX,resid
              !p.multi = 0
              for i=0l,nparams-1l do begin
                 print,'P['+strtrim(i,1)+'] = ',$
                       strtrim(result[i],1),' +/- ',$
                       strtrim(perr[i],1)
              endfor
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
           staramps[j,0] = result[6]
           starampsE[j,0] = perr[6]
           staramps[j,1] = result[7]
           starampsE[j,1] = perr[7]

           acheck = widths[0:j,0]
        endif else begin
           for i=0l,1 do begin
              ;; Go through both sources as long as they're
              ;; sufficiently far apart
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
                 expr = 'voigt_approx(P[0],X - P[1],P[2]) * P[3]'+$
                        '+ eval_legendre(X,P[4:6])'
                 start = [0.1E,posit[i],3E,max(b) * 3E,median(b),0E,0E]
                 nparams = n_elements(start)
                 pi = replicate({fixed:0, limited:[0,0], limits:[0.0E,0.0E]},nparams)
                 pi[0].limited = [1,0]
                 pi[0].limits = [0E,0E] ;; don't let the Voigt parameter be less than 0
                 ;; don't let the star be farther from the slit than the slit length
                 pi[1].limits = [1,1]
                 pi[1].limits = [startx,endx]
                 pi[2].limited = [1,1]
                 pi[2].limits = [-HW/2E,HW/2E] ;; don't let the width be too extreme
                 
                 ;; Estimate outside points for error +/- the baseline
                 ;;                                       half width
                 outp = where((indices LT startx and indices GE startx - HW) or $
                              (indices GE endx   and indices LT endx + HW))
                 yerr = robust_sigma(b[outp]) + fltarr(endx-startx+1)

                 result = mpfitexpr(expr,subarrayX,subarrayY,yerr,start,parinfo=pi,perr=perr,/quiet)
                 yfit = expression_eval(expr,subarrayX,result)
                 resid = yfit - b[startx:endx]

;                 result = gaussfit(subarrayX,subarrayY,A,nterms=5)
;                 Widths[j,i] = A[2]
;                 starLocations[j,i] = A[1]

                 ;;save the results
                 Widths[j,i] = result[2]
                 WidthsE[j,i] = perr[2]
                 starlocations[j,i] = result[1]
                 starlocationsE[j,i] = perr[1]
                 voigts[j,i] = result[0]
                 voigtsE[j,i] = perr[0]
                 staramps[j,i] = result[3]
                 starampsE[j,i] = perr[3]
              endelse
              if keyword_set(showplot) then begin
                 !p.multi = [0,1,2]
                 
                 if i EQ 0 then begin
                    zsubarrayX = subarrayX
                    zresult = result
                    zyfit = yfit
                    zresid = resid
                 endif else begin
                    plot,indices,b,xtitle='Y position (px)',$
                         ytitle='Flux'
                    oplot,subarrayX,yfit,color=mycol('blue'),thick=thickline
                    oplot,zsubarrayX,zyfit,color=mycol('blue'),thick=thickline
                    
                    al_legend,['Measured','Model'],linestyle=[0,0],$
                              color=[!p.color,mycol('blue')],$
                              thick=[thinline,thickline]
                    al_legend,[cgGreek('sigma')+' = '+$
                               string(result[2],format='(F6.2)')+' +/- '+$
                               string(perr[2],format='(F6.2)'),$
                               'a = '+$
                               string(result[0],format='(F6.2)')+' +/- '+$
                               string(perr[0],format=('(F6.2)'))],$
                              /right
;                         'a = '+string(result[
                    plot,indices,b,xtitle='Y position (px)',$
                         ytitle='Residual',yrange=threshold(resid),/nodata
                    oplot,subarrayX,resid
                    oplot,zsubarrayX,zresid
                    for i=0l,nparams-1l do begin
                       print,'P['+strtrim(i,1)+'] = ',$
                             strtrim(result[i],1),' +/- ',$
                             strtrim(perr[i],1)
                    endfor
                    !p.multi = 0

;                    oplot,subarrayX,result,color=mycol('blue'),thick=thickline
;     oplot,subarrayX,A[0] * exp(-((subarrayX - A[1])/A[2])^2) + A[3]
;                                               + A[4] *
;                                               SubarrayX,thick=2,color=mycol('blue')
                    ;; stop
;                    oplot,zsubarrayX,zresult,color=mycol('lblue'),thick=thickline
;                    if j GE 50 then stop
                 endelse
              endif
           endfor
           
        endelse
        
     endfor
     if not keyword_set(specific) then begin
        save,Widths,starLocations,voigts,$
             widthsE,starlocationsE,voigtsE,$
             starAmps,$
             filename='data/state_parameters/widths/widths_'+$
             specfileListNamePrefix+'.sav'
     endif
  endif

  restore,'data/state_parameters/widths/widths_'+$
          specfileListNamePrefix+'.sav'
  save,Widths,starLocations,voigts,$
       widthsE,starlocationsE,voigtsE,$
       starAmps,$
       filename='data/prof_widths.sav'

  if keyword_set(psplot) then begin
     device, /close
;     cgPS2PDF,plotprenm+'.eps'
;     spawn,'convert -density 300% '+plotprenm+'.pdf'+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
     !p.thick=1
     !x.thick=1
     !y.thick=1
  endif

end


