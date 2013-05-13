pro ev_hypertrain,parinfo=pi,findsig=findsig
;; Finds the Gaussian Process hyper-parameters for a training data set
;; uses the tnmin program to optimize the hyper-parameters

  readcol,'data/simulated_series/simser.txt',xt,yt,$
          skipline=1,format='(F,F)'

  ;; Try evaluating the hypertrain parameters
  pstart = [0.4,100,1]
;  pi =
;  replicate({value:0.D,fixed:0,limited:[1,0],limits:[0.D,0.D]},3)

;;; Try a simple task of fitting data to a line to test code
  xt = indgen(100)
  yt = indgen(100) + randomn(0,100) * 5
  stop
  
  pi = replicate({value:0.D,fixed:0,limited:[0,0],limits:[0.D,0.D]},3)
  pi(2).fixed = 1 ;; fix the white noise sigma parameter for now
  pi(*).value = pstart

  FUNCTARGS = {x:xt,y:yt} ;; these will be ased
  a1 = ev_leval(pstart,x=xt,y=yt)

  fitp = tnmin('ev_leval',FUNCTARGS=FUNCTARGS,parinfo=pi,/autoderivative)

;  return,fitparams
end
