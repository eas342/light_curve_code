pro train_on_resids
;; fits the hyperparameters to a set of residuals

  readcol,'data/cleaned_tim_ser/timeser_1.43um_.txt',$
;  readcol,'data/cleaned_tim_ser/timeser_0.91um_.txt',$
;  readcol,'data/cleaned_tim_ser/timeser_1.08um_.txt',$
          phase,fl,flerr,modelfl,resid
;  stop
  phaseFact = max(phase)
;  phaseFact = 1D
  phase = phase / phaseFact
  multfac = 1E
;  residFact = 0.1 * max(resid)
  residFact = 1D
  resid = resid / residFact * multfac

  custpi = replicate({value:0.D,fixed:0,limited:[1,0],limits:[1D-5,0.D],$
                      step:0.001*multfac},3)

;  custstart=[1,3,1]

  custstart = [0.01 * multfac,0.05,0.27 * multfac]
;  custstart = [0.01,0.05,0.27]

   ParamArr = ev_hypertrain(phase,resid,pstart=custstart,$
                           parinfo=custpi)
;   stop
   ;; Re-scale the variables
   finalParamArr = [ParamArr[0] * residFact,$
                    ParamArr[1] * phaseFact,$
                    ParamArr[2] * residFact]


   save,finalParamArr,$
     filename='data/fitted_hyperparams.sav'

end
