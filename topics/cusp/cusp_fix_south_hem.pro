
; use 1 time, fix sign error of ionr for southern hemisphere.
; New cusp_read_plot_data already generates the correct outputs.
pro cusp_fix_south_hem, rootdir

    fns = file_search(rootdir+'/*_all_data.tplot', count = nfn)
    if nfn eq 0 then stop
    
    for i = 0, nfn-1 do begin
;        print, fns[i]
;        store_data, '*', /delete
;        tplot_restore, filename = fns[i]
;        get_data, 'event_info', tmp, dat
;        id = dat.id
;        get_data, 'scidat', tmp, dat
;        hem = dat.polar.cusp.entry.ilat
;        if hem gt 0 then continue
;        print, id
;        dat.ionr *= -1
;        tplot_save, '*', filename = fns[i]
    endfor
end

rootdir = shomedir()+'/Google Drive/works/data/cusp'
;cusp_fix_south_hem, rootdir
end