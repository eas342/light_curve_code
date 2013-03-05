pro run_multiple_bins,differential=differential,quadfit=quadfit
;;differential - shows the differential time series
  night=['dec23','jan04']
  suffix=['png','pdf']

  for i=0,n_elements(night)-1 do begin
     for j=2,13 do begin ;; cycle through # of wavelength bins
        case night[i] of
           'dec23': compile_spec,/optimal,nwavbins=j,/dec23
           'jan04': compile_spec,/optimal,nwavbins=j
        endcase
        plot_tim_ser,/fitcurve,/freelimb,/quadfit,/offtranserr
        if keyword_set(differential) then begin
           for k=0,1 do begin
              plot_tim_ser,/fitcurve,/freelimb,quadfit=quadfit,/offtranserr,$
                           timebin=25,/differential,fixrad=k,/pngcopy,/psplot
              for l=0,n_elements(suffix)-1 do begin
                 if k EQ 1 then comment='_fixrad' else comment=''
                 spawn,'mv plots/spec_t_series/tser_1.22.'+suffix[l]+$
                       ' plots/spec_t_series/tser_1.22_02pt_'+night[i]+comment+'.'+suffix[l]
                 spawn,'mv plots/spec_t_series/tser_2.01.'+suffix[l]+$
                       ' plots/spec_t_series/tser_2.01_02pt_'+night[i]+comment+'.'+suffix[l]
              endfor
           endfor
        endif else begin
           plot_rad_vs_wavl,/psplot
           spawn,'mv plots/rad_vs_wavl.png plots/rad_vs_wavl_'+night[i]+'_'+strtrim(j,1)+'pts.png'
           spawn,'mv plots/rad_vs_wavl.pdf plots/rad_vs_wavl_'+night[i]+'_'+strtrim(j,1)+'pts.pdf'
        endelse
     endfor
  endfor
end
