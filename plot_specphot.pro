pro plot_specphot,divbymodel=divbymodel,usebin=usebin,removelin=removelin,$
                  psplot=psplot,individual=individual,skipInitialize=skipInitialize,$
                  timebin=timebin,backg=backg,custYmargin=custYmargin,$
                  differential=differential,filter=filter,noNorm=noNorm,$
                  backratio=backratio
;; Makes an image of the spectrophotometry to get a visual sense of
;; the transit
;; divbymodel -- divide the image by the nominal transit model
;; usebin -- use the wavelength bins
;; individual -- specifies which of the inidividual spectra to look at
;; backg -- use the background instead of the star flux
;; backratio -- look at the ratio fo the background between the two stars
;; removelin -- remove the linear trend in each time series
;; psplot -- makes a postscript plot
;; skipInitialize -- skips running plot_tim_ser to run faster
;; timebin -- uses time-binned data
;; custYmargin -- for use by double_specphot for shrinking Y margin
;; differential -- use a differential lightcurve instead of absolute
;; filter -- use a digital filter to take out the broad spectral shapes
;; noNorm -- don't normalize the spectrum

  ;; get the time info

  if not keyword_set(skipInitialize) then begin
     plot_tim_ser
  endif
  restore,'data/timedata.sav'
  restore,'data/specdata.sav'

  ntime = n_elements(tplot)

  case 1 of
     n_elements(individual) NE 0: begin
        nwavs = n_elements(lamgrid)
        xydivspec = transpose(flgrid[*,individual-1,*],[0,2,1])
        wavrange = [lamgrid[0],lamgrid[nwavs-1l]]
     end
     keyword_set(usebin): begin
        nwavs = n_elements(bingrid)
        xydivspec = binfl
        wavrange = [bingrid[0],bingrid[nwavs-1]]
     end
     keyword_set(backratio): begin
        nwavs = n_elements(lamgrid)
        xydivspec = transpose(backdiv[*,0,*],[0,2,1])
        wavrange = [lamgrid[0],lamgrid[nwavs-1l]]
     end
     keyword_set(backg): begin
        nwavs = n_elements(lamgrid)
        xydivspec = transpose(backgrid[*,backg-1,*],[0,2,1])
        wavrange = [lamgrid[0],lamgrid[nwavs-1l]]
     end
     keyword_set(snrRatio): begin
        nwavs = n_elements(lamgrid)
        xydivspec = transpose((flgrid[*,0,*]/errgrid[*,0,*])/(flgrid[*,1,*]/errgrid[*,1,*]),[0,2,1])
        wavrange = [lamgrid[0],lamgrid[nwavs-1l]]
     end
     else: begin
        nwavs = n_elements(lamgrid)
        xydivspec = transpose(divspec[*,0,*],[0,2,1])
        wavrange = [lamgrid[0],lamgrid[nwavs-1l]]
     end
  endcase
  ;; Make a median spectrum to divide out
  meddivspec = fltarr(nwavs)
  for i=0l,nwavs-1l do begin
     meddivspec[i] = median(xydivspec[i,*])
  endfor
  replicatedspec = rebin(meddivspec,nwavs,ntime)
  if keyword_set(individual) OR keyword_set(noNorm) then begin
     xypic = xydivspec
  endif else begin
     xypic = xydivspec / replicatedspec
     ;; Normalize by median spectrum
  endelse

  if keyword_set(removelin) then begin
     ;;throw away all n_sigma events before de-trending
     firstCutSig = 12E

     offp = where(tplot LT hstart OR tplot GT hend)
     for i=0l,nwavs-1l do begin

        if total(finite(xypic[i,offp])) EQ 0 then goodp = [-1] else begin
           rstdoff = robust_sigma(xypic[i,offp])
           medoff = median(xypic[i,offp])
           
           goodp = where(abs(xypic[i,*] - medoff) LE firstCutSig * rstdoff and $
                         (tplot LT hstart OR tplot GT hend),complement=throwaways)
        endelse
        if goodp NE [-1] then begin
;           if throwaways NE [-1] then begin
;              tclip1 = tplot[throwaways]
;              yclip1 = y[throwaways]
;           endif
;           yfull = y
           ytemp = xypic[i,goodp]
;           divbycurveclip1 = divbycurve[goodp]
           yerr = xypic[i,goodp]
           tplottemp = tplot[goodp]
;           airmass = airmass[goodp]
           ;; fit result to a robust line
           rlinefit = robust_linefit(tplottemp,ytemp)
           xypic[i,*] = xypic[i,*] / (rlinefit[0] + tplot * rlinefit[1])
        endif
;           ;; divide by the line to flatten out
;           yflat = divbycurveclip1 / yfit
     endfor
  endif

  if keyword_set(filter) then begin
     xypic = convol(xypic,digital_filter(0.15,0.3,50,10))
  endif

  if keyword_set(divbymodel) then begin
     ;; divide all time series by the transit model
     ymodel = quadlc(tplot,planetdat.p,planetdat.b_impact,$
                     u1parm,u2parm,planetdat.a_o_rstar)
     replicatedmodel = rebin(ymodel,ntime,nwavs)
     rebinmodel = transpose(replicatedmodel,[1,0])
     xypic = xypic / rebinmodel
  endif
  if keyword_set(differential) then begin
     restore,'data/cleaned_curve.sav'
     replicatedtrend = rebin(mainCurve,ntime,nwavs)
     rebinTrend = transpose(replicatedTrend,[1,0])
     xypic = xypic / rebinTrend
  endif

  if keyword_set(divbymodel) then begin
     ColorRange = [0.995E,1.005E]
;  endif else ColorRange = [0.98E,1.005E]
  endif else begin
     readcol,'param_input/specphot_range.txt',lowRange,highRange,$
             skipline=1,format='(F,F)'
     ColorRange = [lowRange[0],highRange[0]]
  endelse 

  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotnmpre = 'plots/specphot_images/specphot_image'
     device,encapsulated=1, /helvetica,$
            filename=plotnmpre+'.eps',bits_per_pixel=8
     device,xsize=10, ysize=9,decomposed=1,/color
  endif

  loadct,1
  if not keyword_set(psplot) and not keyword_set(skipInitialize) then window,0
  if n_elements(custYmargin) EQ 0 then custYmargin=[4,2]
  plotimage,xypic,range=ColorRange,$
            imgxrange=wavrange,$
            imgyrange=[tplot[0],tplot[ntime-1]],$
            xtitle='Wavelength (um)',$
            ytitle='Orbital Phase',$
            charsize=1,$
            xmargin=[9,12],ymargin=custYmargin
  ;; Back an outline for thick axes
  axis,xaxis=1,xrange=!x.crange,color=mycol('black'),xstyle=1,xthick=4,xtickname=replicate(' ',7)
  axis,xaxis=1,xrange=!x.crange,color=mycol('orange'),xstyle=1,xtickname=replicate(' ',7)
  axis,xaxis=0,xrange=!x.crange,color=mycol('black'),xstyle=1,xthick=4,xtickname=replicate(' ',7)
  axis,xaxis=0,xrange=!x.crange,color=mycol('orange'),xstyle=1,xtickname=replicate(' ',7)
  axis,yaxis=1,yrange=!y.crange,color=mycol('black'),ystyle=1,ythick=4,ytickname=replicate(' ',7)
  axis,yaxis=1,yrange=!y.crange,color=mycol('orange'),ystyle=1,ytickname=replicate(' ',7)
  axis,yaxis=0,yrange=!y.crange,color=mycol('black'),ystyle=1,ythick=4,ytickname=replicate(' ',7)
  axis,yaxis=0,yrange=!y.crange,color=mycol('orange'),ystyle=1,ytickname=replicate(' ',7)

  loadct,0
  ;; Show ingress and egress
  plots,[wavrange[0],wavrange[1]],[hstart,hstart],color=mycol('black'),linestyle=2,thick=6
  plots,[wavrange[0],wavrange[1]],[hstart,hstart],color=mycol('yellow'),linestyle=2,thick=3
  plots,[wavrange[0],wavrange[1]],[hend,hend],color=mycol('black'),linestyle=2,thick=6
  plots,[wavrange[0],wavrange[1]],[hend,hend],color=mycol('yellow'),linestyle=2,thick=3
  loadct,1

  ;; Choose a set of parameters to pass on to the fits file for the header
  
  keepFitsParams = ['TELESCOP','INSTRUME','OBSERVER','OBJECT','COMMENT',$
                    'DATE_OBS','ITIME','NDR','SLIT','GRAT','GFLT','AFOC',$
                    'PLATE_SC','RA','DEC','EPOCH','EXPTIME']
  nkeepParams = n_elements(keepFitsParams)

  ;; Additional parameters
  addedParams = ['XSTART','XFINISH','DELTAX',$
                 'YSTART','YFINISH','DELTAY']
  addedValues = [lamgrid[0],lamgrid[nwavs-1l],lamgrid[1]-lamgrid[0],$
                 tplot[0],tplot[ntime-1l],tplot[1]-tplot[0]]

  addedComments = [' / Microns start wavelength',$
                   ' / Microns End wavelength',$
                   ' / Wavelength step (microns)',$
                  ' / Orbital Phase start',$
                  ' / Orbital Phase end',$
                  ' / Orbital Phase step']

  nadded = n_elements(addedParams)
  outHeader = strarr(nkeepParams + nadded)
  for i=0l,nkeepParams-1l do begin
     HeaderInd = strpos(header,keepFitsParams[i])
     GoodInd = where(headerInd EQ 0)
     outHeader[i] = Header[goodInd[0]]
  endfor
  for i=0l,nadded-1l do begin
     outHeader[i + nKeepParams] = string(addedParams[i],format='(A-8)') + '=' +$
                                  string(addedValues[i],format='(F21.7)') + $
                                  string(addedComments[i],format='(A-50)')
  endfor

  ;; Save as a FITS image
  fitsNamePre = 'data/specphot'
  if keyword_set(divbymodel) then fitsNamePre = fitsNamePre+'_divided'
  if keyword_set(usebin)     then fitsNamePre = fitsNamePre+'_bin'
  if keyword_set(removelin)  then fitsNamePre = fitsNamePre+'_lin_detrend'
  
  ;; Generate a primary FITS header
  mkhdr,simpleHeader,xypic
  nsimple = n_elements(simpleHeader)
  fullHeader = [simpleHeader[0l:nsimple-4l],outHeader,$
                string('END',format='(A-60)'),string(' ',format='(A-60)')]

  writefits,fitsNamePre+'.fits',xypic,fullHeader

  ;; make an image for the legend
  Ncolors =  256l
  legrow = findgen(Ncolors)*(ColorRange[1]-ColorRange[0])/float(Ncolors)+ColorRange[0]
  legimg = transpose(rebin(legrow,Ncolors,3))
  plotimage,legimg,imgyrange=ColorRange,$
            range=ColorRange,$
            ytitle='Relative Flux',$
            xstyle=4,/noerase,position=[0.95,!y.window[0],0.98,!y.window[1]]
  if keyword_set(psplot) then begin
     device,/close
     cgPS2PDF,plotnmpre+'.eps'
     spawn,'convert -density 450% '+plotnmpre+'.pdf '+plotnmpre+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

end
