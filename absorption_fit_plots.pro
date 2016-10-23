pro absorption_fit_plots
    ;; Plots the telluric absorption line fits
    
    restore,'data/specdata.sav'
    ;; array with shifts-3D [feature #, aperture #, image #]
    restore,'data/telluric/tellshifts.sav'
    readcol,'data/telluric/fit_features.txt',$
            twave,twidth,skipline=1
    
    nap = n_elements(apkey);; number of apertures
    nimg= n_elements(utgrid) ;; number of images
    nfeatures = n_elements(twave)
    
    imgnumarr = findgen(nimg) ;; image number array
    
    legendNames1D = strarr(nfeatures * nap)
    legendNames2D = reform(legendNames1D,nfeatures,nap)

    colarr1D = myarraycol(nfeatures)
    colarr2D = rebin(colarr1D,nfeatures,nap)
    sizeArr1D = findgen(nap) + 1.0
    sizeArr2D = rebin(transpose(sizeArr1D),nfeatures,nap)
    for i=0l,nap-1l do begin
       for j=0l,nfeatures-1l do begin
          ;;for j=0l,0l do begin
          yplot = tellshiftpx[j,i,*]
          if (i EQ 0) AND (j EQ 0) then begin
             plot,imgnumarr,yplot,/nodata,$
                  xtitle='Image #',ytitle='Shift (px)',$
                  yrange=[-5,5]
          endif
          oplot,imgnumarr,yplot,color=colarr2D[j,i],psym=4,symsize=sizeArr2D[j,i]
          legendNames2D[j,i] = 'Nap= '+strtrim(i,1)+' Feature='+strtrim(twave[j])
       endfor 
    endfor
    al_legend,legendNames2D,color=colarr2D,psym=4,symsize=sizeArr2D
    
 end
