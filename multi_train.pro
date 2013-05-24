pro multi_train
;; generates several sets of data, and fits the hyperparameters

Ntrys = 100l
Nparams = 3l

ParamArr = fltarr(Nparams,Ntrys) ;; array for the fitted parameters

ParamTrue = [10E,20E,1E]

for i=0l,Ntrys-1l do begin
   simulate_series,nrealizations=i+1,theta=[ParamTrue[0],ParamTrue[1]],$
                   sigma=ParamTrue[2],Npoints=100l
   readcol,'data/simulated_series/simser.txt',x,y,$
           skipline=1,format='(D,D)'

   ParamArr[*,i] = ev_hypertrain(x,y)
endfor

save,ParamArr,ParamTrue,$
     filename='data/training_series.sav'

end
