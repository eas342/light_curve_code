pro twotick_labels,xtickvals,ytickvals,$
                   xorient=xorient,yorient=yorient,$
                   ymid=ymid,xmid=xmid,$
                   noX=noX,noY=noY,$
                   ytitle=ytitle,xtitle=xtitle
;; Takes the original tick names and only shows the first and last
;; tick labels to avoid over-crowding the plot
;; [XY]tickvals are the tick values found from plot,[XY]tick_get
;; [XY]orient - is the optional keyword to angle the x axis text
;; [XY]mid - does the middle two tick labels (not the first and last)
;; no[XY] - skips a set of axis labels
;; [XY]title - puts in a ytitle left of the words
     numxticks = n_elements(xtickvals)
     numyticks = n_elements(ytickvals)
     dataperYpix = (!y.crange[1] - !y.crange[0]) /((!y.window[1] - !y.window[0]) * !d.y_vsize)
     dataperXpix = (!x.crange[1] - !x.crange[0]) /((!x.window[1] - !x.window[0]) * !d.x_vsize)

     ;; Do the X axis
     if not keyword_set(noX) then begin
        if n_elements(xorient) EQ 0 then xorient=0
        xleft = xtickvals[0]
        xright = xtickvals[numxticks-1l]
        xbottom = !y.crange[0]- 1.5 * dataperYpix * !D.Y_CH_SIZE
        
        xyouts,xleft, xbottom,string(xtickvals[0],format='(G0)'),alignment=0.5,orientation=xorient
        xyouts,xright,xbottom,string(xtickvals[numxticks-1l],format='(G0)'),alignment=0.5,orientation=xorient
     endif

     if not keyword_set(noY) then begin
        ;; Do the Y axis
        if n_elements(yorient) EQ 0 then begin
           yorient=0 ;; default is tilted text centered on value
           yalign = 1
        endif else yalign = 1
        if keyword_set(ymid) then pt2 = floor(numyticks/2E) else pt2 = numyticks-1l
        yleft = !x.crange[0] - 0.3 * dataperXpix * !D.Y_CH_SIZE
        ytop = ytickvals[0] - 0.5 * dataperYpix * !D.Y_CH_SIZE
        
        ybottom = ytickvals[pt2] - 0.5 * dataperYpix * !D.Y_CH_SIZE

        xyouts,yleft,ytop,   string(ytickvals[0],format='(G0)'),$
               alignment=yalign,orientation=yorient
        xyouts,yleft,ybottom,string(ytickvals[pt2],format='(G0)'),$
               alignment=yalign,orientation=yorient
     endif

     if n_elements(ytitle) NE 0 then begin ;; Y title
        yleft = !x.crange[0] - 4.5 * dataperXpix * !D.Y_CH_SIZE
        ycent = (!y.crange[0] + !y.crange[1])/2E
        xyouts,yleft,ycent,ytitle,orientation=90,alignment=0.5
     endif

end
