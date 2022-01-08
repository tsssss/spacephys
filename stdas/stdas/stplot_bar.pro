pro stplot_bar, var, val, position = pos
    
    if n_elements(pos) ne 4 then $
        tplot, var[0], get_plot_position = pos
    get_data, var[0], tmp, tmp, alimit = lims
    yrange = lims.yrange
    plot, [0,1], yrange, position = pos, /normal, $
        _extra = lims, /noerase, /nodata, xtickformat = '(A1)'
    oplot, [0,1], [val,val]
end