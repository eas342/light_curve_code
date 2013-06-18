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
     
     ;; If one of the numbers is close to 10^-15 and should be rounded
     ;; then round
     roundp = where(abs(xtickvals) LT 1E-16,complement=noround)
     if n_elements(roundp) EQ 1 and roundp NE [-1] then begin
        ;; there should either be exactly 1 point to round if the
        ;; numbers go like -0.1,1E-16,0.1,0.2
        xtickvals[roundp] = float(round(xtickvals[roundp]*1000l)/1000E)
     endif
     roundp = where(abs(ytickvals) LT 1E-16,complement=noround)
     if n_elements(roundp) EQ 1 and roundp NE [-1] then begin
        ;; there should either be exactly 1 point to round if the
        ;; numbers go like -0.1,1E-16,0.1,0.2
        ytickvals[roundp] = float(round(ytickvals[roundp]*1000l)/1000E)
     endif

     ;; Do the X axis
     if not keyword_set(noX) then begin
        if n_elements(xorient) EQ 0 then begin
           xorient=0
           xalign=0.5
        endif else xalign=1
        if keyword_set(xmid) then pt2 = floor(numxticks/2E) else pt2 = numxticks-1l
        xleft = xtickvals[0]
        xright = xtickvals[pt2]
        xbottom = !y.crange[0]- 1.2 * dataperYpix * !D.Y_CH_SIZE
        
        xyouts,xleft, xbottom,string(xtickvals[0],format='(G0)'),$
               alignment=xalign,orientation=xorient
        xyouts,xright,xbottom,string(xtickvals[pt2],format='(G0)'),$
               alignment=xalign,orientation=xorient
        ;; Extend the tick marks
;        plots,replicate(xtickvals[0],2),!y.crange[0] + [-1E,1E] * dataperYpix * !D.y_CH_SIZE * 0.8,thick=2
;        plots,replicate(xtickvals[pt2],2),!y.crange[0] + [-1E,1E] * dataperYpix * !D.y_CH_SIZE * 0.8,thick=2
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
        ;; Extend the tick marks
;        plots,!x.crange[0] + [-0.3E,1E] * dataperXpix * !D.y_CH_SIZE * 0.8,thick=2,replicate(ytickvals[0],2)
;        plots,!x.crange[0] + [-0.3E,1E] * dataperXpix * !D.y_CH_SIZE * 0.8,thick=2,replicate(ytickvals[pt2],2)
     endif

     if n_elements(ytitle) NE 0 then begin ;; Y title
        yleft = !x.crange[0] - $
                (2 + max(strlen(string(ytickvals,format='(G0)')))) * 0.7 * dataperXpix * !D.Y_CH_SIZE
        ycent = (!y.crange[0] + !y.crange[1])/2E
        xyouts,yleft,ycent,ytitle,orientation=90,alignment=0.5
     endif

     if n_elements(xtitle) NE 0 then begin ;; X title
        xcent = (!x.crange[0] + !x.crange[1])/2E
        if n_elements(xorient) EQ 0 then lowerfact=3 else begin
;           lowerfact= (2 +
;           max(strlen(string(xtickvals,format='(G0)')))) * 0.7
           lowerfact = 6.3
        endelse
        xbottom = !y.crange[0] - lowerfact * dataperYpix * !D.Y_CH_SIZE
        xyouts,xcent,xbottom,xtitle,alignment=0.5
     endif

end
