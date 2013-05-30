pro twotick_labels,xtickvals,ytickvals,xorient=xorient
;; Takes the original tick names and only shows the first and last
;; tick labels to avoid over-crowding the plot
;; [XY]tickvals are the tick values found from plot,[XY]tick_get
;; xorient - is the optional keyword to angle the x axis text
     numxticks = n_elements(xtickvals)
     numyticks = n_elements(ytickvals)
     dataperYpix = (!y.crange[1] - !y.crange[0]) /((!y.window[1] - !y.window[0]) * !d.y_vsize)
     dataperXpix = (!x.crange[1] - !x.crange[0]) /((!x.window[1] - !x.window[0]) * !d.x_vsize)

     ;; Do the X axis
     if n_elements(xorient) EQ 0 then xorient=0
     xleft = xtickvals[0]
     xright = xtickvals[numxticks-1l]
     xbottom = !y.crange[0]- 1.5 * dataperYpix * !D.Y_CH_SIZE

     xyouts,xleft, xbottom,string(xtickvals[0],format='(G0)'),alignment=0.5,orientation=xorient
     xyouts,xright,xbottom,string(xtickvals[numxticks-1l],format='(G0)'),alignment=0.5,orientation=xorient

     ;; Do the Y axis
     yleft = !x.crange[0] - 0.3 * dataperXpix * !D.Y_CH_SIZE
     ytop = ytickvals[0]
     ybottom = ytickvals[numyticks-1l]
     
     xyouts,yleft,ytop,   string(ytickvals[0],format='(G0)'),alignment=0.5,orientation=90
     xyouts,yleft,ybottom,string(ytickvals[numyticks-1l],format='(G0)'),alignment=0.5,orientation=90

end
