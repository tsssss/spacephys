;+
; PROCEDURE: pstplot.pro
;
; PURPOSE: Wrapper for TPLOT which sends to tplot to the printer or a 
;            postscript file rather than the screen
;
; INPUTS:  DATANAMES  -> Names or IDs of tplot quantities (same as TPLOT
;                          first argument). (Optional)
;
; KEYWORDS: Same as TPLOT keywords, plus ...
;           FILENAME   -> File name to save output to
;           NO_COLOR   /  Plots in greyscale
;           FULL       /  Plots in full page format rather than the proportions
;                         of the default IDL window
;           CHSIZE     -> Base Character size to use
;           TBARS      -> array of times for time bars
;           BARCOLORS  -> colors to use for time bars
;           BARPANELS  -> panels to put time bars in, if not whole plot
;                          (use -1 in array for whole plot)
;           NOCLOSE    /  Don't close plot
;           CLOSE      /  Used after /NOCLOSE to finish plot
;           [XY]SIZE   -> Size of plot in inches
;           [XY]OFFSET -> Offset of plot in inches
;
; CALLING SEQUENCE: pstplot,/full
;                   pstplot,filename='out.ps',/full
;
; NOTES: Scaling is designed to work off of default IDL window.  If you 
;         rescale the IDL window used in TPLOT, PSTPLOT may not scale properly.
;        The easiest method is to use regular TPLOT to setup the plot as
;         you want to see it, and then just do PSTPLOT(,/FULL).  
;     ** This is a wrapper for TPLOT, so if you use any other TPLOT keywords
;         they will effect the plot on the screen as well as the PS output.
;        If a file name is not given, output will go to the default printer.
;
; CREATED BY: John P. Dombeck 1/24/2001
;
; MODIFICATION HISTORY:
;
;  Initial Version - 1/24/2001 by John P. Dombeck
;  07/25/01-J. Dombeck       Added CHSIZE keyword
;  11/02/01-J. Dombeck       Added TBARS,BARCOLORS,NOCLOSE,CLOSE keywords
;  02/15/02-J. Heisserer     Added [XY]SIZE,[XY]OFFSET,Changed to inches
;  04/09/02-J. Dombeck       Added BARPANELS keyword
;  06/26/02-J. Dombeck       Auto send to default printer if filename not given
;  09/16/03-J. Dombeck       Auto append ".ps" to filename if not given
;-
;INCLUDED MODULES:
;   pstplot
;
;LIBRARIES USED:
;   tplot_com
;
;DEPENDANCIES
;   tplot
;
;-



;*** MAIN *** : * PSTPLOT *

pro pstplot,datanames,filename=filename,no_color=no_color,full=full, $
            chsize=chsize,tbars=tbars,barcolors=barcolors,barpanels=barpanels,$
            noclose=noclose,close=close,xsize=xsize,ysize=ysize,$
            xoffset=xoffset,yoffset=yoffset,$
            axis2 = axis2, _EXTRA=ex

common pstplot_static,old_charsize

@tplot_com.pro

; Close plot

  if keyword_set(close) then begin
    device,/close
    case !version.os_family of
        'unix': set_plot, 'x'
        'Windows': set_plot, 'win'
    endcase
    !p.charsize=old_charsize
    if n_elements(filename) eq 0 then spawn,'lpr idl.ps'
    return
  endif

; Check tbar matching

  nbars=n_elements(tbars)
  nclrs=n_elements(barcolors)
  if nbars ne nclrs and nclrs ne 0 then begin
    message,'TBARS / BARCOLORS mis-match',/cont
    return
  endif
  npnls=n_elements(barpanels)
  if nbars ne npnls and npnls ne 0 then begin
    message,'TBARS / BARPANLES mis-match',/cont
    return
  endif


; Save old value of chsize

  old_charsize=!p.charsize
  if not keyword_set(chsize) then chsize=old_charsize
  if chsize eq 0. then chsize=1.
 

; Set device

  set_plot,'PS'
  device,filename=filename,/encapsulate
  
  if not keyword_set(no_color) then device,/color,bits=32  $
    else device,color=0,bits=32
  

; Format output

  if keyword_set(full) then begin
    !p.charsize=1.2*chsize
    device,ysize=10.2,yoffset=0.4,xsize=7.7,xoffset=0.4,/inches
  endif else begin
    !p.charsize=0.7*chsize
    device,ysize=9.0,yoffset=1.00,xsize=6.5,xoffset=1.0,/inches
  endelse
  if keyword_set(xsize) then device,xsize=xsize,/inches
  if keyword_set(ysize) then device,ysize=ysize,/inches
  if keyword_set(xoffset) then device,xoffset=xoffset,/inches
  if keyword_set(yoffset) then device,yoffset=yoffset,/inches

; Plot

  tplot,datanames,_EXTRA=ex
  if n_elements(axis2) ne 0 then $
      tplot2ndaxis, axis2.var, axis2.range, title = axis2.title


; Time bars

  if nbars gt 0 then begin
    if npnls gt 0 then begin
      if nclrs eq 0 then begin
        for xx=0,nbars-1 do begin
          if barpanels[xx] eq -1 then timebar,tbars[xx] $
                                 else timebar,tbars[xx],window=barpanels[xx]
        endfor
      endif else begin
        for xx=0,nbars-1 do begin
          if barpanels[xx] eq -1 then timebar,tbars[xx],color=barcolors[xx] $
                                 else timebar,tbars[xx],color=barcolors[xx],$
                                        window=barpanels[xx]
        endfor
      endelse
    endif else begin
      for xx=0,nbars-1 do begin
        if nclrs eq 0 then timebar,tbars[xx] $
                      else timebar,tbars[xx],color=barcolors[xx]
      endfor
    endelse
  endif


; Clean up

  if not keyword_set(noclose) then pstplot,/close,filename=filename

return
end       ;*** MAIN *** : * PSTPLOT *

