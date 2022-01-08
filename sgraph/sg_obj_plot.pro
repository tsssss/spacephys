function sg_og_pos2view, pos, xr, yr
    rec = dblarr(4)
    rec[2] = (xr[1]-xr[0])/(pos[2]-pos[0])
    rec[0] = xr[0]-rec[2]*pos[0]
    rec[3] = (yr[1]-yr[0])/(pos[3]-pos[1])
    rec[1] = yr[0]-rec[3]*pos[1]
    return, rec
end

; a scene is a window, view is a plot. see mg_ogscene_example

pro sg_obj_plot, x0, y0, $
    aspect_ratio = aspect, $
    background_color = bg_color, $
    bakcground_transparency = bg_trans, $
    clip = clip, $
    color = color, $
    eye = eyedist, $
    histogram = histogram, $
    linestyle = linestyle, $
    max_value = maxval, $
    min_value = minval, $
    position = pos, $
    rgb_table = rgbct, $
    stairstep = stairstep, $
    symbol = symbol, $
    thick = thick, $
    title = title, $
    transparency = trans, $
    window = win, $
    xrange = xr, $
    yrange = yr, $
    xsize = xsize, $
    ysize = ysize
    
    ; [xy]range and position.
    if n_elements(xr) eq 0 then xr = (n_params() eq 1)? $
        [0,n_elements(x0)]: [min(x0,/nan),max(x0,/nan)]
    if n_elements(yr) eq 0 then yr = (n_params() eq 1)? $
        [min(x0,/nan),max(x0,/nan)]: [min(y0,/nan),max(y0,/nan)]
    if n_elements(pos) eq 0 then pos = [0.13d,0.13,0.92,0.88]
    viewpos = sg_og_pos2view(pos, xr, yr)

    ; [xy]size and window.
    if n_elements(xsize) eq 0 then xsize = !d.x_size
    if n_elements(ysize) eq 0 then ysize = !d.y_size
    if n_elements(win) eq 0 then begin
        container = obj_new('IDL_Container')
        view = obj_new('IDLgrView', viewplane_rect = viewpos)
        model = obj_new('IDLgrModel')
        win = obj_new('IDLgrWindow', retain = 2, dimensions = [xsize,ysize])
    endif
    font = obj_new('IDLgrFont', 'times')
    
    plot = obj_new('IDLgrPlot', x0, y0)
    tmp = {exact:1,notext:1}
    axb = sg_og_axis(0, location = [xr[0],yr[0]], range = xr)
    axt = sg_og_axis(0, location = [xr[0],yr[1]], range = xr, tickdir = 1, /notext)
    ayl = sg_og_axis(1, location = [xr[0],yr[0]], range = yr)
    ayr = sg_og_axis(1, location = [xr[1],yr[0]], range = yr, tickdir = 1, /notext)
    
    view->add, model
    model->add, plot
    model->add, axb
    model->add, axt
    model->add, ayl
    model->add, ayr
    
    container->add, win
    container->add, view
    container->add, font
    
    win->draw, view
    
    stop
    obj_destroy, container
end

x = findgen(201)
y = sin(2*!dpi/25*x)*exp(-0.02*x)
sg_obj_plot, x, y
end
