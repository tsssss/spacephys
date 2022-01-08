
pro cg_zbuffer_brain

    ; initial settings.
    dev0 = 'x'
    xsize = 500
    ysize = 500
    set_plot, 'z'
    !p.font = 1
    device, decomposed = 0, set_pixel_depth = 24
    device, set_resolution = [xsize,ysize]
    help, /device
    
    ; prepare data.
    fn = filepath(subdir = ['examples','data'], 'head.dat')
    head = bytarr(80,100,57)
    openr, lun, fn, /get_lun
    readu, lun, head
    free_lun, lun
    
    xpt = 60 & ypt = 64 &  zpt = 20
    ximg = reform(head[xpt,*,*])
    yimg = reform(head[*,ypt,*])
    zimg = reform(head[*,*,zpt])
    s = size(head, /dimensions)
    xs = s[0]-1 & ys = s[1]-1 & zs = s[2]-1
    xplane = [[xpt,0,0],[xpt,0,zs],[xpt,ys,zs],[xpt,ys,0],[xpt,0,0]]
    yplane = [[0,ypt,0],[0,ypt,zs],[xs,ypt,zs],[xs,ypt,0],[0,ypt,0]]
    zplane = [[0,0,zpt],[xs,0,zpt],[xs,ys,zpt],[0,ys,zpt],[0,0,zpt]]
    scale3, xrange = [0,xs], yrange = [0,ys], zrange = [0,zs]
    
    ncolors = 233
    black = '000000'x
    white = 'ffffff'x
    
    erase, white
    loadct, 0, ncolors = ncolors
    polyfill, xplane[*,0:3], /t3d, pattern = ximg, $
        image_coord = [[0,0],[0,zs],[ys,zs],[ys,0]], /image_interp
    loadct, 1, ncolors = ncolors
    polyfill, yplane[*,0:3], /t3d, pattern = yimg, $
        image_coord = [[0,0],[0,zs],[xs,zs],[xs,0]], /image_interp
    loadct, 3, ncolors = ncolors
    polyfill, zplane[*,0:3], /t3d, pattern = zimg, $
        image_coord = [[0,0],[xs,0],[xs,ys],[0,ys]], /image_interp
    
    thick = 2
    plots, xplane, /t3d, thick = thick, color = black
    plots, yplane, /t3d, thick = thick, color = black
    plots, zplane, /t3d, thick = thick, color = black

    tvlct, r, g, b, /get
    snap = tvrd()
    snap24 = tvrd(true = 3)
    help, snap, snap24
    set_plot, dev0
    window, 0, xsize = xsize, ysize = ysize
    tvimage, snap, /nointerp
    window, 1, xsize = xsize, ysize = ysize
    tvimage, snap24, /nointerp
    write_png, 'brain.png', snap, r, g, b
    write_jpeg, 'brain24.jpg', snap24, true = 3
end