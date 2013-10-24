pro prepare_kepler_data,psplot=psplot
;; Takes a phase folded Kepler light curve and prepares it for a
;; functional fit
;; psplot - save a plot of the orignal folded data

  ;; Get average light curve
  readcol,'data/phase_folded_kepler_all_kic1255.txt',kbinnum,kphaseS,kfluxS,format='(F,F,F)',/silent

  ;; Change phase to 0
  x = kphaseS - 1.0

  phaseX = x
  save,phaseX,kfluxS,filename='data/kepler_curves/phase_folded_kepler_all_kic1255.sav'

  ;; Make the Out of transit points flat & normalize
  
  y = kfluxS 

  outp = where((x GT -0.6 and x LE -0.1) or x GT 0.3)
  avgout = mean(y[outp])
  
  ;; re-normalize
  y = y/avgout

  ;; set the plot
  if keyword_set(psplot) then begin
     set_plot,'ps'
     !p.font=0
     plotprenm = 'plots/kep_curves/kep_avg'
     device,encapsulated=1, /helvetica,$
            filename=plotprenm+'.eps'
           device,xsize=5, ysize=5,decomposed=1,/color
  endif
  ;; Do a normal plot first
  plot,kphaseS-1.0,y,ystyle=16,xrange=[-0.5,0.5],$
       xtitle='Orbital Phase',xstyle=1,psym=10,charsize=0.4
  if keyword_set(psplot) then begin
     device, /close
     cgPS2PDF,plotprenm+'.eps'
     spawn,'convert -density 400% '+plotprenm+'.pdf '+plotprenm+'.png'
     device,decomposed=0
     set_plot,'x'
     !p.font=-1
  endif

  ;; make all out of transit points = 1
  y[outp] = 1.0E

  ;; truncate points below -0.5 and above 0.5
  truncp = where(x GE -0.5 and x LE 0.5)
  x = x[truncp]
  y = y[truncp]
  plot,x,y,ystyle=16

;  z = (x - 0.05)/0.01
;  z2 = (x - 0.07)/0.01
;  f = -exp(-z^2)*0.005 -exp(-z2^2)*0.005 + 1
  
  x2 = x - 0.05
;  f = -5.38 * gaussian(29E * x2 * gaussian(2.91 * x2,[1,0,1]) -
;  0.251,[1,0,1])
;  f = 0.0409242E - $
;      3.32628E *atan(0.000934522E, atan(x2)^2 - atan(0.00554009E *x2,1.23707E + 38.4843E *x2) - 0.00618158E * atan(x2))

;  f = -3.373E * atan(0.001195E + 0.01908E *x2, x2 *tan(2.368E *x2))
  f = -3.2 * atan(0.00093,x2^2 - x2 * atan(0.012,1.1 + 26E * x2))
;  f = -3.778E *atan(0.000929E, x2^2 - atan(-0.0002431E, 1.224E + 37.88E *x2))
  f = f * 0.001E + 1E
  oplot,x,f,color=mycol('yellow')

end
