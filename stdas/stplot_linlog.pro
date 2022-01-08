;+
; ytitlepos, in ycharsize, default is 5.
; set lastpanel to show the var_labels.
;-
pro stplot_linlog, var0, position = tpos, labels = labs, $
    linyrange = linyr, linytickv = linytkv, linyminor = linymnr, $
    logyrange = logyr, logytickv = logytkv, logyminor = logymnr, $
    ypans = ypans, ytitlepos = dytit, $
    lastpanel = lastp, $
    _extra = ex


    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size


    if n_elements(tpos) eq 0 then tpos = sgcalcpos()
    if n_elements(linymnr) eq 0 then linymnr = 5
    if n_elements(logymnr) eq 0 then logymnr = 5
    if n_elements(dytit) eq 0 then dytit = 4
    if n_elements(linestyle) eq 0 then linestyle = 1
    
    ; get initial settings.
    get_data, var0, uts, dat, limits = lim

    
    ; divide up position.
    x0 = tpos[0]
    x2 = tpos[2]
    y0 = tpos[1]
    y2 = tpos[3]
    if n_elements(ypans) ne 2 then ypans = [1,1]
    ratio0 = double(ypans[1])/total(ypans)
    y1 = y0+(y2-y0)*ratio0
    
    pos1 = [x0,y1,x2,y2]
    pos2 = [x0,y0,x2,y1]
    
    
    tvar = var0+'1'
    store_data, tvar, uts, dat, limits = lim
    options, tvar, 'xstyle', 1+4
    options, tvar, 'yrange', logyr
    options, tvar, 'ylog', 1
    options, tvar, 'yticks', 2
    options, tvar, 'ystyle', 1
    options, tvar, 'xtickformat', '(A1)'
    options, tvar, 'ytitle', ''
    options, tvar, 'labels', ''
    if n_elements(logytkv) ne 0 then begin
        options, tvar, 'ytickv', logytkv
        options, tvar, 'yticks', n_elements(logytkv)-1
    endif
    options, tvar, 'yminor', logymnr
    tplot, tvar, position = pos1, /noerase, /novtitle, _extra = ex
    plots, pos1[[0,2]], pos1[1]+[0,0], /normal, linestyle = linestyle
    plot, [0,1],[0,1], position = pos1, /nodata, /noerase, xstyle = 5, ystyle = 5
    axis, xaxis = 1, xtickformat = '(A1)', xticklen = 0, xticks = 1, xminor = 0, xstyle = 1


    tvar = var0+'2'
    store_data, tvar, uts, dat, limits = lim
    options, tvar, 'xstyle', 1+8
    options, tvar, 'yrange', linyr
    options, tvar, 'ylog', 0
    options, tvar, 'ystyle', 1
    options, tvar, 'ytitle', ''
    options, tvar, 'labels', ''
    if n_elements(linytkv) ne 0 then begin
        options, tvar, 'ytickv', linytkv
        options, tvar, 'yticks', n_elements(linytkv)-1
    endif
    options, tvar, 'yminor', linymnr
    
    if keyword_set(lastp) then begin
        novtit = 0
        nottk = 0
    endif else begin
        novtit = 1
        nottk = 1
        options, tvar, 'xtickformat', '(A1)'
    endelse
    
    tplot, tvar, position = pos2, /noerase, novtitle = novtit, $
        nouttick = nouttk, _extra = ex

    tx = x0-ychsz*dytit
    ty = (y0+y2)*0.5
    xyouts, tx, ty, /normal, alignment = 0.5, orientation = 90, lim.ytitle
    
    labs = lim.labels
    cols = lim.colors
    nlab = n_elements(labs)
    tx = x2+xchsz*1
    if nlab eq 1 then ys = (y0+y2)*0.5-ychsz*0.5 else begin
        ys = smkarthm(y2-ychsz, y0, (y0-y2)/nlab, 'dx')
    endelse
    ys = ys+ychsz*0.2
    for i = 0, nlab-1 do $
        xyouts, tx, ys[i], /normal, alignment = 0, labs[i], color = cols[i]


    ; restore initial settings.
    store_data, var0+['1','2'], /delete



end

erase
device, decomposed = 0
loadct2, 43

pre0 = 'rbspb_'
tvar = pre0+'pf_fac_mat'
get_data, tvar, uts, dat, limit = lims

tvar = tvar+'_map'
store_data, tvar, uts, dat*500, limit = lims
options, tvar, 'labels', ['b','e','n']
stplot_linlog, tvar, position = sgcalcpos(), ytitlepos = 2, $
    linyr = [-5,5], logyr = [5,100], logytickv = [10,100], linytickv = [-5,0,5]
end
