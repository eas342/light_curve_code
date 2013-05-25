pro train_on_resids
;; fits the hyperparameters to a set of residuals

  readcol,'data/cleaned_tim_ser/timeser_1.43um_.txt',$
;  readcol,'data/cleaned_tim_ser/timeser_0.91um_.txt',$
          phase,fl,flerr,modelfl,resid

  phaseFact = max(phase)
  phase = phase / phaseFact
  residFact = 0.1 * max(resid)
  resid = resid / residFact
  
;   ParamArr = ev_hypertrain(phase,resid,pstart=[0.1,0.001 *
;   1E4,0.05])

  custpi = replicate({value:0.D,fixed:0,limited:[1,0],limits:[1D-8,0.D],$
                      step:0.05},3)

   ParamArr = ev_hypertrain(phase,resid,pstart=[1,3,1],custStepsize=0.5E,$
                           parinfo=custpi)
;   stop
   ;; Re-scale the variables
   finalParamArr = [ParamArr[0] * residFact,$
                    ParamArr[1] * phaseFact,$
                    ParamArr[2] * residFact]


   save,finalParamArr,$
     filename='data/fitted_hyperparams.sav'

end
