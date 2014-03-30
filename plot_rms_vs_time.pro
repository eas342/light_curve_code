pro plot_rms_vs_time
;; Plot the variation as a function of time

readcol,'data/cleaned_tim_ser/timeser_1.61um_.txt',$
        phase,fl,flRMS,model,resid,format='(F,F,F,F,F)'


plot,phase,flRMS,xrange=[-0.01,0.21],xstyle=1

end
