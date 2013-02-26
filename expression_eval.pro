function expression_eval,expr,X,P
  ;; Evaluates the kind of expression that is put into MPFIT
  junk = execute('y= '+expr)
  return,y
end
