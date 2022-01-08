;+
; Use triangulation to get a 2-D movie of B. Doesn't work well because the distribution of magnetometers is bad.
;-

;---Get data.
    seed = -1l
    x = randomu(seed, 41)
    y = randomu(seed, 41)
    z = exp(-3*((x-0.5)^2+(y-0.5)^2))

    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1
    mlat_range = einfo.mag_mlat_range
    mlon_range = einfo.mag_mlon_range
    mag_mlons = einfo.mag_mlons
    mag_mlats = einfo.mag_mlats

    mag_var = 'thg_dbh'
    get_data, mag_var, times, dbh, sites, limit=lims
    x = mag_mlons
    y = mag_mlats
    
;---Prepare the plot.
    sgopen, 0, xsize=8, ysize=4
    plot, x, y, /nodata, xrange=mlon_range, yrange=mlat_range, xstyle=1, ystyle=1
    plots, x, y, psym=1, color=sgcolor('red')


;---Create the trianglues from xy.
    triangulate, x, y, triangles
    s = size(triangles, /dimensions)
    ntriangle = s[1]
    
    for ii=0, ntriangle-1 do begin
        thistriangle = [triangles[*,ii],triangles[0,ii]]
        plots, x[thistriangle], y[thistriangle], color=sgcolor('teal')
    endfor
    
    ct = 71
    nlevel = 51
    levels = smkarthm(-800,800,nlevel,'n')
    colors = smkarthm(0,250,nlevel,'n')
    for ii=0, nlevel-1 do colors[ii] = sgcolor(colors[ii],ct=ct)

    for ii=0, n_elements(times)-1, 12 do begin
        
        ; Create the grid.
        gridded_data = trigrid(x,y,reform(dbh[ii,*]), triangles, $
            xgrid=xvector, ygrid=yvector)
    
        contour, gridded_data, xvector, yvector, xstyle=1, ystyle=1, $
            xrange=mlon_range, yrange=mlat_range, position=sgcalcpos(1), $
            nlevel=nlevel, levels=levels, /fill, c_colors=colors

        xyouts, 0.1,0.1, /normal, time_string(times[ii]), color=sgcolor('black')
        wait, 0.02
        
    endfor
end
