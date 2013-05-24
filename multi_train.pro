pro multi_train
;; generates several sets of data, and fits the hyperparameters

Ntrys = 100l
Nparams = 3l

ParamArr = fltarr(Nparams,Ntrys) ;; array for the fitted parameters

thetaTrue = [10E,20E]

for i=0l,Ntrys-1l do begin
   simulate_series,nrealizations=i+1,theta=thetaTrue,Npoints=100l
   readcol,'data/simulated_series/simser.txt',x,y,$
           skipline=1,format='(D,D)'

   ParamArr[*,i] = ev_hypertrain(x,y)
endfor

save,ParamArr,thetaTrue,$
     filename='data/training_series.sav'

end
