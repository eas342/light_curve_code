pro ev_hypertrain,parinfo=pi,findsig=findsig
;; Finds the Gaussian Process hyper-parameters for a training data set
;; uses the tnmin program to optimize the hyper-parameters

  readcol,'data/simulated_series/simser.txt',xt,yt,$
          skipline=1,format='(D,D)'

  ;; Try evaluating the likelihood function
;  pstart = [0.4,100,1]
;  pstart = [0.05,240,1]
  pstart = [10,20,1]

  pi =  replicate({value:0.D,fixed:0,limited:[1,0],limits:[0.D,0.D],$
                 step:0.05},3)

  pi(2).fixed = 1 ;; fix the white noise sigma parameter for now
  pi(*).value = pstart

  FUNCTARGS = {x:xt,y:yt} ;; these will be ased
  a1 = ev_leval(pstart,x=xt,y=yt)

  fitp = tnmin('ev_leval',FUNCTARGS=FUNCTARGS,parinfo=pi,/autoderivative)
;  stop
;  return,fitparams
end
