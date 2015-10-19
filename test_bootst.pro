pro test_bootst,seed1=seed1,showp=showp
;; Generates some data to test the bootstrap method

expr = 'poly(X,P)'
if n_elements(seed1) EQ 0 then seed1 = 0

np = 128
xtest = findgen(np)/float(np -1l) * 5E
true = [0,0,1E]

ytest = poly(xtest,true) + randomn(seed1,np)
yerr = fltarr(np) + 1E

dat = struct_arrays(create_struct('XTEST',xtest,'YEST',ytest))



start = [0,0,0]
;result =
;mpfitexpr(expr,xtest,ytest,yerr,start,parinfo=pi,perr=perr,/quiet)
result = boot_mp(expr,xtest,ytest,yerr,start,parinfo=pi,perr=perr,mperr=mperr)
model = expression_eval(expr,xtest,result)

ev_oplot,dat,xtest,model,gparam=gparam

sigFormat = '(A8,3F10.3)'
print,'Input: '
print,'True',true,format=sigFormat
print,'InSig',fltarr(n_elements(true)),format=sigFormat
print,' '

;print,'Uncertainties:'
;print,'MPFIT',result,format=sigFormat
;print,'MP+/-',mperr,format=sigFormat
;print,'Boot+/-',perr,format=sigFormat
;

if keyword_set(showp) then genplot,dat,gparam=gparam


end
