pro avg_radii,totsets=totsets,statistics=statistics,psplot=psplot,$
              median=median
;;totsets -- optional keyword to specify the number of total sets of
;;           data to avg, default is 2
;; statistics -- look at the statistics of the variations
;; psplot - make postscript plot
;; median - find a median instead of a weighted average

  if n_elements(totsets) EQ 0 then totsets=2

  ;; read in the first radius versus wavelength file
  radfile = choose_file(searchDir='radius_vs_wavelength',$
                           filetype='.txt')
  ;; Search the radius_vs_wavlength directory for .txt files

  readcol,radfile,wavl,wavlsize,rad,rade,skipline=1,format='(F,F,F)'

  nwavs = n_elements(wavl)
  avgwavl = wavl
  avgwavlsize = wavlsize
  totrad = fltarr(nwavs,totsets)
  totrade = fltarr(nwavs,totsets)
  avgrad = fltarr(nwavs)
  avgrade = fltarr(nwavs)

  totrad[*,0] = rad
  totrade[*,0] = rade
     
  for i=1l,totsets-1l do begin
     print,'Choose Additional file ',strtrim(i-1l,1)
     file2 = choose_file(searchDir='radius_vs_wavelength',filetype='.txt')
     readcol,file2,wavl2,wavl2size,rad2,rade2,skipline=1,format='(F,F,F)'
     if n_elements(wavl2) NE nwavs then begin
        print,'Different number of wavelength bins in files!'
        return
     endif else begin
        assert,total(wavl2-avgwavl),'=',0,'Different Wavelength!'
        assert,total(wavl2size-avgwavlsize),'=',0,'Different Wavelength Sizes!'
        totrad[*,i] = rad2
        totrade[*,i] = rade2
     endelse
  endfor

  ;; Do a weighted average
  for j=0l,nwavs-1l do begin
     weights = 1E / totrade[j,*]^2
     avgrad[j] = total(weights * totrad[j,*]) / total(weights)
     avgrade[j] = 1E / sqrt(total(weights))
  endfor

  ;; Do a median, if asked
  if keyword_set(median) then begin
     avgRad = median(totrad,dimension=2)
  endif

  forprint,avgwavl,avgwavlsize,avgrad,avgrade,$
           textout='radius_vs_wavelength/avg_rp_rs.txt',$
           comment='# Wavelength(um)  Bin size   Rp/R*    Rp/R* Error'

  if keyword_set(statistics) then begin
     if keyword_set(psplot) then begin
        set_plot,'ps'
        !p.font=0
        plotprenm = 'plots/radius_distribution/rad_distrib'
        device,encapsulated=1, /helvetica,$
               filename=plotprenm+'.eps'
        device,xsize=14, ysize=10,decomposed=1,/color
     endif


     ;; Look at the statistics of the data points
     fullavg = rebin(avgrad,nwavs,totsets)
     diffavg = (totrad - fullavg)/totrade
     binsize=0.5
     yhist = histogram(diffavg,locations=xhist,binsize=binsize)
     plot,xhist,yhist,psym=10,$
          xtitle='(R!Dp!N/R!D*!N - Mean)/Error',$
          ytitle='Count'
     ;; Finish histogram bars
     fp = n_elements(yhist)-1l
     oplot,[xhist[0] - binsize/2E,xhist[0] - binsize/2E,xhist[0]],$
           [!y.crange[0],yhist[0],yhist[0]]
     oplot,[xhist[fp],xhist[fp] + binsize/2E,xhist[fp] + binsize/2E],$
           [yhist[fp],yhist[fp],!y.crange[0]]


     ;; Gaussian
     xgauss = (indgen(256)/256E - 0.5E) * 10E
     dx = xgauss[1] - xgauss[0]
     ygauss = exp(- xgauss^2/2E)
     ygauss = ygauss / total(ygauss * dx) * total(yhist * binsize)
     oplot,xgauss,ygauss,color=mycol('blue'),linestyle=2

     al_legend,['Fitted Radii','Gaussian PDF'],$
            linestyle=[0,2],color=[!p.color,mycol('blue')]

     if keyword_set(psplot) then begin
        device, /close
        cgPS2PDF,plotprenm+'.eps'
        spawn,'convert -density 350% '+plotprenm+'.pdf '+plotprenm+'.png'
        device,decomposed=0
        set_plot,'x'
        !p.font=-1
     endif

  endif

end  
  
