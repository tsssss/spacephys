;+
; Type: procedure.
; Purpose: Advanced plot wrapper.
; Parameters:
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   2014-04-04, Sheng Tian, create.
;   2015-11-05, Sheng Tian, use sg_[prep,free]_data, use struct as plot input.
;-

pro sglogplot, x0, y0, $
    position = position, $
    xlim = xlim0, ylim = ylim0, xrange = xrange, yrange = yrange, $
    xtitle = xtitle, ytitle = ytitle, title = title, psym = psym, $
    noxtick = noxtick, noytick = noytick, _extra = extra

    ; convert input to pointer.
    sg_prep_data, x0, y0, px = px, py = py, nparam = nparam

    ; **** initial settings.
    opt = {$
        color:sgcolor('black'), $
        background:sgcolor('white'), $
        noerase:1, $
        name: 'options'}
    if keyword_set(noxtick) then opt = create_struct('xtickformat','(A1)', opt)
    if keyword_set(noytick) then opt = create_struct('ytickformat','(A1)', opt)

    if n_elements(position) eq 0 then position = sgcalcpos()
    if n_elements(psym) eq 0 then psym = 0

    ; determine range.
    if not keyword_set(xrange) then xrange = [-1e2,1e3]
    if not keyword_set(yrange) then yrange = [-1e2,1e3]
    
    ; if [xy]lim are within [xy]range.
    if n_elements(xlim0) eq 0 then xlim = 1d1 else xlim = double(xlim0)
    if n_elements(ylim0) eq 0 then ylim = 1d1 else ylim = double(ylim0)
    if n_elements(xlim) eq 1 then xlim*= [-1,1]
    if n_elements(ylim) eq 1 then ylim*= [-1,1]
        
    xneg = where(*px lt 0, nxneg)
    yneg = where(*py lt 0, nyneg)
    
    sglogplot_calcpos, position[0], position[2], xrange[0], xrange[1], $
        xlim[0], xlim[1], x10, x11, x12, x13
    sglogplot_calcpos, position[1], position[3], yrange[0], yrange[1], $
        ylim[0], ylim[1], y10, y11, y12, y13
    
    ; plot box.
    plot, [x10,x13],[y10,y13], /nodata, $
        position = position, $
        ticklen = 0, xtickformat = '(A1)', ytickformat = '(A1)', $
        xstyle = 1, ystyle = 1, _extra = opt
    plots, [x10,x13], [y11,y11], linestyle = 1, psym = 0, _extra = opt
    plots, [x10,x13], [y12,y12], linestyle = 1, psym = 0, _extra = opt
    plots, [x11,x11], [y10,y13], linestyle = 1, psym = 0, _extra = opt
    plots, [x12,x12], [y10,y13], linestyle = 1, psym = 0, _extra = opt
    
    ; plot data.
    ; x>0, y>0.
    idx = where(*px gt xlim[1] and *py gt ylim[1], cnt)
    txr = [xlim[1],xrange[1]]
    tyr = [ylim[1],yrange[1]]
    if cnt gt 0 then begin
        plot, (*px)[idx], (*py)[idx], $
            position = [x12,y12,x13,y13], psym = psym, $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr, tyr, /nodata, $
            position = [x12,y12,x13,y13], $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
    axis, xaxis = 1, xlog = 1, xstyle = 1, xtickformat = '(A1)', _extra = opt
    axis, yaxis = 1, ylog = 1, ystyle = 1, ytickformat = '(A1)', _extra = opt
    
    ; x>0, y<0.
    idx = where(*px gt xlim[1] and *py lt ylim[0], cnt)
    txr = [xlim[1],xrange[1]]
    tyr =-[yrange[0],ylim[0]]
    if cnt gt 0 then begin
        plot, (*px)[idx],-(*py)[idx], $
            position = [x12,y10,x13,y11], psym = psym, $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr, tyr, /nodata, $
            position = [x12,y10,x13,y11], $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
    axis, xaxis = 0, xlog = 1, xstyle = 1, _extra = opt
    axis, yaxis = 1, ylog = 1, ystyle = 1, ytickformat = '(A1)', _extra = opt
    
    ; x<0, y>0.
    idx = where(*px lt xlim[0] and *py gt ylim[1], cnt)
    txr =-[xrange[0],xlim[0]]
    tyr = [ylim[1],yrange[1]]
    if cnt gt 0 then begin
        plot,-(*px)[idx], (*py)[idx], $
            position = [x10,y12,x11,y13], psym = psym, $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr, tyr, /nodata, $
            position = [x10,y12,x11,y13], $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
    axis, xaxis = 1, xlog = 1, xstyle = 1, xtickformat = '(A1)', _extra = opt
    axis, yaxis = 0, ylog = 1, ystyle = 1, _extra = opt
    
    ; x<0, y<0.
    idx = where(*px lt xlim[0] and *py lt ylim[0], cnt)
    txr =-[xrange[0],xlim[0]]
    tyr =-[yrange[0],ylim[0]]
    if cnt gt 0 then begin
        plot,-(*px)[idx],-(*py)[idx], $
            position = [x10,y10,x11,y11], psym = psym, $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr, tyr, /nodata, $
            position = [x10,y10,x11,y11], $
            xlog = 1, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
    axis, xaxis = 0, xlog = 1, xstyle = 1, xtickformat = '("-",I0)', _extra = opt
    axis, yaxis = 0, ylog = 1, ystyle = 1, ytickformat = '("-",I0)', _extra = opt

    ; mid-right.
    idx = where(*px gt xlim[1] and *py ge ylim[0] and *py le ylim[1], cnt)
    txr = [xlim[1],xrange[1]]
    tyr = ylim
    if cnt gt 0 then begin
        plot, (*px)[idx], (*py)[idx], $
            position = [x12,y11,x13,y12], psym = psym, $
            xlog = 1, xrange = txr, $
            ylog = 0, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr, tyr, $
            position = [x12,y11,x13,y12], $
            xlog = 1, xrange = txr, $
            ylog = 0, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
;    axis, yaxis = 1, ylog = 0, yminor = 1, yticks = 2, ystyle = 1, ytickformat = '(A1)', _extra = opt
    
    ; mid-left.
    idx = where(*px le xlim[0] and *py ge ylim[0] and *py le ylim[1], cnt)
    txr =-[xrange[0],xlim[0]]
    tyr = ylim
    if cnt gt 0 then begin
        plot,-(*px)[idx], (*py)[idx], $
            position = [x10,y11,x11,y12], psym = psym, $
            xlog = 1, xrange = txr, $
            ylog = 0, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr, tyr, /nodata, $
            position = [x10,y11,x11,y12], $
            xlog = 1, xrange = txr, $
            ylog = 0, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
;    axis, yaxis = 0, ylog = 0, yminor =1 , yticks = 2, ystyle = 1, ytickformat = '(A1)', _extra = opt
    
    ; up-mid.
    idx = where(*px ge xlim[0] and *px le xlim[1] and *py ge ylim[0], cnt)
    txr = xlim
    tyr = [ylim[1],yrange[1]]
    if cnt gt 0 then begin
        plot, (*px)[idx], (*py)[idx], $
            position = [x11,y12,x12,y13], psym = psym, $
            xlog = 0, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr, tyr, /nodata, $
            position = [x11,y12,x12,y13], $
            xlog = 0, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
;    axis, xaxis = 1, xlog = 0, xminor = 1, xticks = 2, xstyle = 1, xtickformat = '(A1)', _extra = opt
    
    ; down-mid.
    idx = where(*px ge xlim[0] and *px le xlim[1] and *py le ylim[1], cnt)
    txr = xlim
    tyr =-[yrange[0],ylim[0]]
    if cnt gt 0 then begin
        plot, (*px)[idx],-(*py)[idx], $
            position = [x11,y10,x12,y11], psym = psym, $
            xlog = 0, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif else begin
        plot, txr,-tyr, $
            position = [x11,y10,x12,y11], $
            xlog = 0, xrange = txr, $
            ylog = 1, yrange = tyr, $
            xstyle = 5, ystyle = 5, _extra = opt
    endelse
;    axis, xaxis = 0, xlog = 0, xminor = 1, xticks = 2, xstyle = 1, xtickformat = '(A1)', _extra = opt

    ; other.
    idx = where(*py ge ylim[0] and *py le ylim[1] and *px gt xlim[0] and *px le xlim[1], cnt)
    txr = xlim
    tyr = ylim
    if cnt gt 0 then begin
        plot, (*px)[idx], (*py)[idx], $
            position = [x11,y11,x12,y12], psym = psym, $
            xlog = 0, xrange = xlim, $
            ylog = 0, yrange = ylim, $
            xstyle = 5, ystyle = 5, _extra = opt
    endif
        
    
        
    ; **** cleanup and restore original settings.
    sg_free_data, px, py, x0 = x0, y0 = y0, nparam = nparam

end

nrec = 1000
x = (findgen(nrec+1)/nrec-0.5)*2*!const.pi
y = 6.5*sin(x)
erase, color = sgcolor('white')
sglogplot, x, y, xrange = [-10,10], yrange = [-10,10], xlim = .5, ylim = .5, xtitle = 'hehe', ytitle = 'haha'
end
