pro kic1255_mult
;; Multiplies the scaling factor times the KSC (Kepler short cadence)
;; transit depth to display the transit depth

  readcol,'radius_vs_wavelength/radius_vs_wavl.txt',$
          bingridmiddle,binsizes,plrad,plrade,format='(F,F,F,F)',skipline=1,/silent

  readcol,'transit_info/kic1255_ksc_depth.txt',mult,format='(F)',skipline=1,/silent

  depth = plrad * mult[0]
  depthErr = plrade * mult[0]

  depthNm = 'radius_vs_wavelength/depth_vs_wavl.txt'
  ;; save the depth data
  forprint,bingridmiddle,binsizes,depth,depthErr,$
           textout=depthNm,$
           comment='#Wavelength(um) Binsize (um)  Depth(%)  Depth Error',/silent

  spawn,'open -R '+depthNm

end
