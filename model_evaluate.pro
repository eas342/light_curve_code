function model_evaluate,X,P
;; Evaluates a saved model
restore,'data/binned_model.sav'
tabinv,X,wavl,Ieff

return,binModel[Ieff] * P

end
