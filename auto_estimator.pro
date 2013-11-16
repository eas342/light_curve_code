function auto_estimator,C
;; Show the covariance estimator (from William
;; Wei's book - "Time Series Analysis"
  Npoints = n_elements(C[0,*])
  VarZ = total(C)/float(Npoints^2)
  kernArray = transpose(C[0,*])
  karray = findgen(Npoints)
  fNpoints = float(Npoints)
  AutoEstimator = kernArray - karray/fNpoints * kernarray - (fNpoints - karray)* VarZ/fNpoints 
  return,AutoEstimator        

end
