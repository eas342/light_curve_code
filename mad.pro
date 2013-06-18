function mad,X
;; Median Absolute deviation
return,median(abs(X - median(X)))
end
