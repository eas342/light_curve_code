function symbol_colors,n,type=type,bback=bback
  ;; This function easily makes arrays of symbols, lines or colors for
  ;; plotting multiple things
  ;; the type keyword sets what to return (symbol, line, color)
  ;; the bback keyword chooses colors that look good on a black
  ;; background. Otherwise, it chooses colors that are good on a white background
  case type of
     'psym': choices = [1,2,4,5,6,7] ;;plot symbol choices (avoiding a dot)
     'linestyle': choices = [0,2,3,4,5] ;; line style choices (avoiding dotted line)
     'color': begin
        if keyword_set(bback) then begin
           ;; colors that are good on a black background
           colorOpt = mycol(['white','red','green','lblue','magenta','yellow'$
                             ,'blue','orange','gold','brown','purple','pink'])
        endif else choices = mycol(['black','red','green','magenta','brown'$
                    ,'blue','orange','purple','pink']) ;; color choices on white paper
     end
     else: choices = [1,2,4,5,6,7] ;;plot symbol choices (avoiding a dot)
  endcase
  nchoices = n_elements(choices)
  symbolarray = choices[lindgen(n) mod nchoices]
  return,symbolarray

end
