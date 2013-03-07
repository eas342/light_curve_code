pro plot_specphot,divbymodel=divbymodel,usebin=usebin
;; Makes an image of the spectrophotometry to get a visual sense of
;; the transit
;; divbymodel -- divide the image by the nominal transit model
;; usebin -- use the wavelength bins

  ;; get the compiled spectroscopic data
  restore,'data/specdata.sav'

  ;; get the time info
  plot_tim_ser
  restore,'data/timedata.sav'

  ntime = n_elements(tplot)


  if keyword_set(usebin) then begin
     nwavs = n_elements(bingrid)
     xydivspec = binfl
     wavrange = [bingrid[0],bingrid[nwavs-1]]
  endif else begin
     nwavs = n_elements(lamgrid)
     xydivspec = transpose(divspec[*,0,*],[0,2,1])
     wavrange = [lamgrid[0],lamgrid[nwavs-1l]]
  endelse

  ;; Make a median spectrum to divide out
  meddivspec = fltarr(nwavs)
  for i=0l,nwavs-1l do begin
     meddivspec[i] = median(xydivspec[i,*])
  endfor
  replicatedspec = rebin(meddivspec,nwavs,ntime)
  xypic = xydivspec / replicatedspec
;  xypic = xydivspec

  if keyword_set(divbymodel) then begin
     ;; divide all time series by the transit model
     ymodel = quadlc(tplot,planetdat.p,planetdat.b_impact,$
                     u1parm,u2parm,planetdat.a_o_rstar)
     replicatedmodel = rebin(ymodel,ntime,nwavs)
     rebinmodel = transpose(replicatedmodel,[1,0])
     xypic = xypic / rebinmodel
  endif

  ColorRange = [0.95E,1.01E]

  loadct,1
  window,0
  plotimage,xypic,range=ColorRange,$
            imgxrange=wavrange,$
            imgyrange=[tplot[0],tplot[ntime-1]],$
            xtitle='Wavelength (um)',$
            ytitle='Orbital Phase',$
            charsize=2,font=1
  loadct,0
  ;; Show ingress and egress
  plots,[wavrange[0],wavrange[1]],[hstart,hstart],color=mycol('brown'),linestyle=2,thick=2
  plots,[wavrange[0],wavrange[1]],[hend,hend],color=mycol('brown'),linestyle=2,thick=2
  loadct,1
  
  fileNpre = 'plots/specphot_images/specphot_image'
  write_png,fileNpre+'.png',tvrd(true=1)
  window,1,ysize=300
  ;; make an image for the legend
  legrow = findgen(256)*(ColorRange[1]-ColorRange[0])/float(256)+ColorRange[0]
  legimg = rebin(legrow,256,3)
  plotimage,legimg,imgxrange=ColorRange,$
            range=ColorRange,$
            xtitle='Relative Flux',font=1,charsize=2,$
            ystyle=4
  write_png,'plots/specphot_images/image_key.png',tvrd(true=1)
  loadct,0

end
