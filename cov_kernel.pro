function cov_kernel,x,theta0,theta1
;; Covariance kernel used by likelihood evaluation, simulated series
;; and plot_tim_ser

;; create array identical to x argument
f = x
Argument = - abs(x^2) * theta1^2
bigNeg = where(Argument LT -15D,complement=smallNeg)
if bigNeg NE [-1] then f[bigNeg] = 0D
if smallNeg NE [-1] then f[smallNeg] = theta0 * exp(Argument[smallNeg])/theta1^2

return,f

end
