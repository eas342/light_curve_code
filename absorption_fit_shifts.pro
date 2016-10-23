pro absorption_fit_shifts,showp=showp
;; Fits absorption lines from telluric features to examine wavelength
;; shifting in some more detail


  nsig = 3 ;; sigma clipping for outliers
  restore,'data/specdata.sav'
  readcol,'data/telluric/fit_features.txt',$
          twave,twidth,skipline=1

  customfunc = 'Gaussian(X,P)'

  nap = n_elements(apkey);; number of apertures
  nimg= n_elements(utgrid) ;; number of images
  nfeatures = n_elements(twave)

  pxgrid = findgen(n_elements(lamgrid))
  tellshiftpx = fltarr(nfeatures,nap,nimg)

  pxleft = fltarr(nfeatures)
  pxright = fltarr(nfeatures)
  pxcenter = fltarr(nfeatures)
  for i=0l,nfeatures-1l do begin
     ;; Get the left, right and center pixels of each feature
     pxleft[i] = interpol(pxgrid,lamgrid,twave[i] - twidth[i] * 0.5E)
     pxright[i] = interpol(pxgrid,lamgrid,twave[i] + twidth[i] * 0.5E)
     pxcenter[i] = interpol(pxgrid,lamgrid,twave[i])
  endfor


  for k=0l,nap-1l do begin
  ;; each aperture
     for j=0l,nimg-1l do begin
        ;; each image
        for i=0,nfeatures-1l do begin

           fitp = where((pxgrid GT pxleft[i]) and (pxgrid LE pxright[i]))
           xtofit = pxgrid[fitp]
           ytofit = flgrid[fitp,k,j]
           
           start = fltarr(4)
           thresh = threshold(ytofit,low=0.2,high=0.8,mult=1.0)
           start[0] = thresh[0]
           start[1] = pxcenter[i]
           start[2] = (pxright[i] - pxleft[i]) * 0.3E
           start[3] = thresh[1]
           fitval = ev_robust_poly(xtofit,ytofit,1,nsig=nsig,$
                                   customfunc=customfunc,start=start,$
                                   _extra=ex)
           center = fitval[1]
           
           tellshiftpx[i,k,j] = center - pxcenter[i]
;           if keyword_set(showp) then begin
           if keyword_set(showp) and abs(tellshiftpx[i,k,j]) GT 2 then begin
              edat = create_struct('VERTLINES',center)
              ymodel = expression_eval(customfunc,xtofit,fitval)
              dat = struct_arrays(create_struct('Wavelength',pxgrid,$
                                                'Flux',flgrid[*,k,j]))
              ev_oplot,dat,xtofit,ymodel,gparam=gparam
              genplot,dat,edat,gparam=gparam
              if quit_caught() then begin
                 return
              endif
              print,filen[j]
              stop
           endif
        endfor
     endfor
  endfor

  save,tellshiftpx,filename='data/telluric/tellshifts.sav'
  stop

end
