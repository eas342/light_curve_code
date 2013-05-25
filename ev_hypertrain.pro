function ev_hypertrain,xt,yt,parinfo=pi,findsig=findsig,$
                       pstart=pstart,custStepsize=custStepsize,$
;; Finds the Gaussian Process hyper-parameters for a training data set
;; uses the tnmin program to optimize the hyper-parameters
;; xt and yt are the input x and y coordinates
;; parinfo is the parameter info (specified the same as in mpfit)
;; findsig -- finds the sigma value
;; pstart -- start values of the hyperparameters [theta_0,theta_1 and sigma]
;; custStepsize -- customo step size
;; parinfo - custom parameter info

;  readcol,'data/simulated_series/simser.txt',xt,yt,$
;          skipline=1,format='(D,D)'

  ;; Try evaluating the likelihood function
;  pstart = [0.4,100,1]
;  pstart = [0.05,240,1]
  if n_elements(pstart) EQ 0 then pstart = [10,20,1]
  if n_elements(custStepsize) EQ 0 then custStepsize = 0.05

  if n_elements(pi) EQ 0 then begin
     pi =  replicate({value:0.D,fixed:0,limited:[1,0],limits:[0.D,0.D],$
                      step:custStepsize},3)
     
     ;;  pi(2).fixed = 1 ;; fix the white noise sigma parameter for now
  endif
  pi(*).value = pstart

  FUNCTARGS = {x:xt,y:yt} ;; these will be ased
  a1 = ev_leval(pstart,x=xt,y=yt)

  fitp = tnmin('ev_leval',FUNCTARGS=FUNCTARGS,parinfo=pi,/autoderivative)

  return,fitp
end
