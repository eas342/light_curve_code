pro correlation_plots
;; Sets up the correlation plots for a given night

  plotDir = 'plots/spatial_pos_corr/'
  readcol,'param_input/correlation_params.txt',plotPren,PkeysX,PkeysY,$
          format='(A,A,A)',skipline=1

  nplots = n_elements(plotPren)
  ;; Get the state parameter info
  restore,'data/used_date.sav'
  restore,'data/state_parameters/full_parameters/'+specfileListNamePrefix+'.sav'
  dat = struct_arrays(statePstruct)
  gparam = create_struct('PSYM',1,'NOMARGLEG',1,'FILENAME','',$
                        'PKEYS',['PHASE','FLUXRATIO'],'SERIES','FLUXRATIO',$
                        'ROUNDSER',0.025)
  for i=0l,nplots-1l do begin
     if file_exists('ev_local_pparams.sav') then begin
        restore,'ev_local_pparams.sav' ;; get the latest parameters
     endif
     gparam.filename=plotDir+usedate+plotpren[i]
     gparam.pkeys[0] = PkeysX[i]
     gparam.pkeys[1] = PkeysY[i]
     gparam.TITLES[0:1] = [PkeysX[i],PkeysY[i]]
     print,'X = ',pkeysX[i],' Y = ',PkeysY[i]
     genplot,dat,gparam=gparam
;     if quit_caught() then return
  endfor

end
