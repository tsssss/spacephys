;+
; patch to add convection velocity and spacecraft velocity.
; polar only.
;-

pro cusp_test_patch_convection, id, polar = polar, fast = fast

    if keyword_set(fast) then forpolar = 0 else forpolar = 1
    if keyword_set(polar) then forpolar = 1 else forpolar = 0
    
    store_data, '*', /delete
    re = 6378d
    rad = !dpi/180
    deg = 180/!dpi

    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')

; **** read basic info.
    ; log info, search on conjunction list, then on polar list.
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    loginfo = cusp_read_conjun_list(logfile, event = id)
    if size(loginfo,/type) ne 8 then begin
        logfile = rootdir+'/works/cusp/cusp_list_of_polar_2-4Re.log'
        loginfo = cusp_read_conjun_list(logfile, event = id, /nofast)
    endif

    if size(loginfo,/type) ne 8 then begin
        print, 'no id found ...'
        return
    endif

    print, 'processing event: '+id


    ; info, id, pre0, plottr, cusptr.
    info = forpolar? loginfo.polar: loginfo.fast
    pre0 = forpolar? 'po_': 'fa_'

    plottr = info.plot_time
    if plottr[1] lt plottr[0] then plottr[1]+= 86400d

    cusptr = info.cusp_time
    if cusptr[1] lt cusptr[0] then cusptr[1]+= 86400d
    
    
    
    ; load scidat.
    datdir = rootdir+'/data/cusp'
    infofn = datdir+'/'+id+'_all_data.tplot'
    infovar = 'scidat'
    tplot_restore, filename = infofn
    
    get_data, infovar, scidatt0, scidat
    scinfo = forpolar? scidat.polar: scidat.fast
    
    if ~forpolar then begin
        dis = (scinfo.cusp.entry.dis+scinfo.cusp.exit.dis)*0.5
        vsc = (abs(scinfo.cusp.exit.ilat)-abs(scinfo.cusp.entry.ilat))/$
            (cusptr[1]-cusptr[0])*dis*re*rad
        plots, dis, abs(vsc), psym = 1
        return
    endif
    
    
    dis = (scinfo.cusp.entry.dis+scinfo.cusp.exit.dis)*0.5
    mlt = (scinfo.cusp.entry.mlt+scinfo.cusp.exit.mlt)*0.5
    dst = scidat.dst
        
    vcv = dis^1.5*0.625
    vsc = (abs(scinfo.cusp.exit.ilat)-abs(scinfo.cusp.entry.ilat))/$
        (cusptr[1]-cusptr[0])*dis*re*rad

    vperp = abs(vcv+vsc)
    print, vperp, dis
    plots, dis, vperp, psym = 1

end

plot, [1,6], [1,10], /nodata, ylog = 0

ids = cusp_calc_id(cusp_id('south_imf'),cusp_id('polar_south_imf'))
foreach id, ids do cusp_test_patch_convection, id, /polar
ids = cusp_id('south_imf')
foreach id, ids do cusp_test_patch_convection, id, /fast
end
