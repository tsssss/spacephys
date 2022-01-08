;+
; horizontal: set to plot a horizontal color bar.
; zrange=. Set zrange, linear always.
; position=.
; ct=.
;-
pro sgcolorbar, colors, zrange = zr, ztitle = ztitle, position = pos, $
    zcharsize = zcharsize, zticks = zticks, zminor = zminor, ztickv = ztickv, ztickname = ztickn, $
    horizontal = horizontal, log = log, zticklen=zticklen, ztickformat=ztickformat, $
    _extra = ex

    xr = [0,1]
    if n_elements(zr) eq 0 then message, 'no zrange ...'
    
    if n_elements(pos) eq 0 then pos = sgcalcpos()
    
    if n_elements(ztitle) eq 0 then ztitle = ''
    if n_elements(zcharsize) eq 0 then zcharsize = 0.8
    
    if n_elements(colors) eq 0 then colors = indgen(256)
    ncolor = n_elements(colors)

    p1 = convert_coord(pos[0], pos[1], /normal, /to_device)
    p2 = convert_coord(pos[2], pos[3], /normal, /to_device)
    
    ;cb = [1d,1] # findgen(ncolor)
    cb = [1d,1] # colors
    if keyword_set(horizontal) then cb = transpose(cb)

    sgtv, cb, position=pos, /resize, _extra = ex
    
    if keyword_set(horizontal) then begin
        plot, zr, xr, position = pos, /normal, /nodata, /noerase, $
            color = sgcolor('black'), background = sgcolor('white'), $
            xrange = zr, ystyle = 1, xstyle = 1, $
            xlog = log, xticks=1, xminor=0, $
            yminor = 0, yticks = 1, ytickformat = '(A1)', xtickformat = '(A1)'
        axis, xaxis = 1, xtitle = ztitle, xcharsize = zcharsize, $
            xrange = zr, xstyle = 1, xticks = zticks, xminor = zminor, $
            xtickv = ztickv, xtickname = ztickn, xtickformat=ztickformat, xticklen=zticklen, $
            color = sgcolor('black')
    endif else begin
        plot, xr, zr, position = pos, /normal, /nodata, /noerase, $
            color = sgcolor('black'), background = sgcolor('white'), $
            yrange = zr, xstyle = 1, ystyle = 1, $
            ylog = log, $
            xminor = 0, xticks = 1, xtickformat = '(A1)', ytickformat = '(A1)'
        axis, yaxis = 1, ytitle = ztitle, ycharsize = zcharsize, $
            yrange = zr, ystyle = 1, yticks = zticks, yminor = zminor, $
            ytickv = ztickv, ytickname = ztickn, ytickformat=ztickformat, yticklen=zticklen, $
            color = sgcolor('black')
    endelse

end

zz = [1,2,3,4,5]
zr = [1,5]
colors = bytscl(zz, top = 120)
ct = 3
device, decomposed = 1
sgopen, 0
pos = [0.85,0.1,0.9,0.9]
sgcolorbar, zrange = zr, colors, position = pos, ct = 3
pos = [0.1,0.8,0.8,0.9]
sgcolorbar, zrange = zr, colors, position = pos, ct = 3, /horizontal
end