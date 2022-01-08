;+
; update Poynting flux statistics and plots.
;-
pro cusp_update_pflux_stat, id

    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    
    ; directory to save output data and pdfs.
    figdir = rootdir+'/works/cusp/cusp list conjun'
    datdir = rootdir+'/data/cusp'
    
    ifn = datdir+'/'+id+'_all_data.tplot'
    tplot_restore, filename = ifn
    
    get_data, 'scidat', ut0, info
    poss = sgcalcpos(2, ypad = 4)
    poss = sgcalcpos()
    erase
    
    tinfo = info.fast
    xr = [1,1e4]
    
    nfilter = tinfo.nfilter
    filters = reverse(tinfo.filters[0:nfilter])
    pfluxes = tinfo.sb.fs[0:nfilter]
    
    plot, filters, pfluxes, xrange = xr, xlog = 1, ylog = 1, position = poss[*,0], $
        xtitle = 'Period (s)', ytitle = 'Integrated S!D||!N (mW/m)', /noerase
    
    tinfo = info.polar
    xr = [1,1e4]
    
    nfilter = tinfo.nfilter
    filters = reverse(tinfo.filters[0:nfilter])
    pfluxes = tinfo.sb.fs[0:nfilter]
    
;    plot, filters, pfluxes, xrange = xr, xlog = 1, ylog = 1, position = poss[*,0], $
;        xtitle = 'Period (s)', ytitle = 'Integrated S!D||!N (mW/m)', /noerase
    oplot, filters, pfluxes, linestyle = 2
    
    stop

end

id = '1998_1001_02'
id = '1998_0925_05'
store_data, '*', /delete
cusp_update_pflux_stat, id
end