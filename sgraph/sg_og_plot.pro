function sg_og_pos2view, pos, xr, yr
    rec = dblarr(4)
    rec[2] = (xr[1]-xr[0])/(pos[2]-pos[0])
    rec[0] = xr[0]-rec[2]*pos[0]
    rec[3] = (yr[1]-yr[0])/(pos[3]-pos[1])
    rec[1] = yr[0]-rec[3]*pos[1]
    return, rec
end

; a scene is a window, view is a plot. see mg_ogscene_example

function sg_og_plot, x0, y0, $
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
    viewpos = sg_og_pos2view([0d,0,1,1], xr, yr)

    ; creat view, use normal coord.
    ogview = obj_new('IDLgrView', location = pos[0:1], $
        dimensions = pos[2:3]-pos[0:1], units = 3, viewplane_rect = viewpos)
    ogmodel = obj_new('IDLgrModel')
    font = obj_new('IDLgrFont', 'times')
    ogmodel->add, obj_new('IDLgrPlot', x0, y0)
    ogmodel->add, sg_og_axis(0, location = [xr[0],yr[0]], range = xr)
    ogmodel->add, sg_og_axis(0, location = [xr[0],yr[1]], range = xr, tickdir = 1, /notext)
    ogmodel->add, sg_og_axis(1, location = [xr[0],yr[0]], range = yr)
    ogmodel->add, sg_og_axis(1, location = [xr[1],yr[0]], range = yr, tickdir = 1, /notext)
    ogview->add, ogmodel
    return, ogview
end

x = findgen(201)
y = sin(2*!dpi/25*x)*exp(-0.02*x)
ogscene = obj_new('IDLgrScene')
pos = sgcalcpos(2)
pos[1,1] = 0 & pos[1,3] = 0.6
ogscene->add, sg_og_plot(x, y, position = pos[0,*])
ogscene->add, sg_og_plot(x, y, position = pos[1,*])
ogwindow = obj_new('idlgrwindow', dimensions = [600,400], graphics_tree = ogscene)
ogwindow->draw
end