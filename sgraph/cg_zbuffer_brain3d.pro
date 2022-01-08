
pro cg_zbuffer_brain3d
    
    ; prepare data.
    fn = filepath(subdir = ['examples','data'], 'head.dat')
    head = bytarr(80,100,57)
    openr, lun, fn, /get_lun
    readu, lun, head
    free_lun, lun
    
    shade_volume, head, 40, vertices, polygons, /low
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
    isosurf = polyshade(vertices, polygons, /t3d)
    stop
    
    black = '000000'x
    white = 'ffffff'x
    device, decomposed = 0
    
    erase, white
    loadct, 1
    tv, isosurf

    tvlct, r, g, b, /get
    snap = tvrd()
    snap24 = tvrd(true = 3)
    help, snap, snap24
    set_plot, dev0
    window, 0, xsize = xsize, ysize = ysize
    tv, snap, /nointerp
    window, 1, xsize = xsize, ysize = ysize
    tv, snap24, /nointerp
    write_png, 'brain3d.png', snap, r, g, b
    write_jpeg, 'brain3d24.jpg', snap24, true = 3
end