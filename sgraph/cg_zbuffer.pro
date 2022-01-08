
pro cg_zbuffer

    dev0 = 'x'
    xsize = !d.x_size
    ysize = !d.y_size
    set_plot, 'z'
    !p.font = 1
    device, decomposed = 1, set_pixel_depth = 24
    device, set_resolution = [xsize,ysize]
    help, /device
    
    peak = shift(dist(20,16), 10,8)
    peak = exp(-(peak/5.0)^2)
    saddle = shift(peak,6,0)+shift(peak,-6,0)/2
    
    loadct, 1, ncolors = 100, bottom = 1
    loadct, 3, ncolors = 100, bottom = 101
    black = '000000'x
    white = 'ffffff'x
    
    set_shading, values = [1,100]
    shade_surf, peak, zrange = [0.0,1.2], charsize = 2.0, $
        color = black, background = white
    set_shading, values = [101,200]
    shade_surf, saddle, zrange = [0.0,1.2], /noerase, charsize = 2.0, $
        color = black
        
    tvlct, r, g, b, /get
    snap = tvrd()
    snap24 = tvrd(true = 3)
    help, snap, snap24
    set_plot, dev0
    window, 0, xsize = xsize, ysize = ysize
    tvimage, snap, /nointerp
    window, 1, xsize = xsize, ysize = ysize
    tvimage, snap24, /nointerp
    write_png, 'zbuffer.png', snap, r, g, b
    write_jpeg, 'zbuffer24.jpg', snap24, true = 3
end