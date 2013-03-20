pro plot_divide_and_indiv,index
;; Plots the host star, the reference star and the the division for an
;; individual wavelength (index)

  ;; get the compiled spectroscopic data
  restore,'data/specdata.sav'

  ;; get the time info
  restore,'data/timedata.sav'

  div1 = flgrid[index,0,*]/flgrid[index,1,*]
  plot,tplot,div1/median(div1),yrange=[0.8,1.2],$
       ytitle='Divided Flux',xtitle='Orbital Phase',$
       title=strtrim(lamgrid[index],1)+'um Flux',charsize=2,font=1,$
       ystyle=8,psym=3
  plot,tplot,flgrid[index,0,*]/median(flgrid[index,0,*]),yrange=[0.8,1.5],$
       ystyle=4,charsize=2,xstyle=4,color=mycol('red'),/noerase,psym=3
  plot,tplot,flgrid[index,1,*]/median(flgrid[index,1,*]) + 0.4E,yrange=!y.crange,$
       ystyle=4,charsize=2,xstyle=4,color=mycol('green'),/noerase,psym=3
  axis,yaxis=1,color=mycol('red'),yrange=!y.crange,$
       ytitle='Individual Star Flux/Median + Offset',font=1,charsize=2
  legend,['Host','Reference'],color=mycol(['red','green']),linestyle=[0,0]

end
