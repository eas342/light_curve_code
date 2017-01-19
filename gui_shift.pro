PRO gui_shift_event, ev

WIDGET_CONTROL, ev.ID, GET_UVALUE=uval ;; retrieve button stuff
widget_control, ev.top, get_uvalue= data ;; retrieve the data structure

CASE uval of
    'NEXT': data.index = wrap_mod((data.index + 1l),data.nfile)
    'PREV': data.index = wrap_mod((data.index - 1l),data.nfile)
    'UP': data.interval = data.interval * 3E
    'DOWN': data.interval = data.interval / 3E
    'LEFT': data.npix[data.index] = data.npix[data.index] - data.interval
    'RIGHT': data.npix[data.index] = data.npix[data.index] + data.interval
    'DONE': begin
       restore,'data/used_date.sav'
       forprint,data.filen,data.npix,textout='data/shift_data/manual/manual_shift'+specfileListNamePrefix+'.txt',$
                format='(A,F)',comment='#Filename                      Shift'

       WIDGET_CONTROL, ev.TOP, /DESTROY
       es_cmd_focus
       return
    end
ENDCASE

manual_shift,data.npix[data.index],nindex=data.index,interval=data.interval

widget_control, ev.top, set_uvalue = data ;; save the data structure

END

pro gui_shift
;; Use the keyboard to shift spectra manually

  base = WIDGET_BASE(/ROW) ;; base to store groups of buttons

  ;; Sets up the control buttons
  menuW = widget_button(base,value = 'File',/menu)
  donebutton = WIDGET_BUTTON(menuW, VALUE='Done', UVALUE='DONE',accelerator='Ctrl+D')
  prevImg = WIDGET_BUTTON(menuW, VALUE='Prev', UVALUE='PREV',accelerator='Ctrl+K')
  nextImg = WIDGET_BUTTON(menuW, VALUE='Next', UVALUE='NEXT',accelerator='Ctrl+L')
  lefW = WIDGET_BUTTON(menuW, VALUE='Left', UVALUE='LEFT',accelerator='Left')
  rightW = WIDGET_BUTTON(menuW, VALUE='Right', UVALUE='RIGHT',accelerator='Right')
  upW = WIDGET_BUTTON(menuW, VALUE='Up', UVALUE='UP',accelerator='Up')
  downW = WIDGET_BUTTON(menuW, VALUE='Down', UVALUE='DOWN',accelerator='Down')
  saveW = WIDGET_BUTTON(menuW, VALUE='Save', UVALUE='SAVE',accelerator='Ctrl+S')

  ;; Sets up the data structure
  restore,'data/specdata.sav'
  nfile = n_elements(utgrid)

  ;; Check for a previous saved manual shift array
  restore,'data/used_date.sav'
  previousName = 'data/shift_data/manual/manual_shift'+specfileListNamePrefix+'.txt'
  if file_exists(previousName) then begin
     readcol,previousName,fileName,initialShift,format='(A,F)'
  endif else begin
     initialShift = fltarr(nfile) -1.7E
  endelse

  data = create_struct('NFILE',nfile,'INDEX',0,'INTERVAL',1E,$
                       'NPIX',initialShift,'FILEN',filen)
  
  manual_shift,data.npix[data.index],nindex=data.index,interval=data.interval

  WIDGET_CONTROL, base, /REALIZE
  widget_control, base, set_uvalue = data

  XMANAGER, 'gui_shift', base

end
