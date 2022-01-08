
pro cusp_update_system, ids, load = load, logxls = logxls, movefile = movefile, combine = combine, $
    test = test, default = default, fullupdate = fullupdate

    syncdatdir = shomedir()+'/Google Drive/works/data/cusp'
    syncfigdir = shomedir()+'/Google Drive/works/works/cusp/cusp list conjun'

    locdatdir = shomedir()+'/cusp/data'
    locfigdir = shomedir()+'/cusp'


    ; default is to update data and figures in sync data dir and sync fig dir.
    datdir = syncdatdir
    figdir = syncfigdir

    ; group the sync fig dir into local fig dir.
    idir = syncfigdir
    odir = locfigdir

    ; log for excel form.
    xlsfn = locfigdir+'/cusp_excel.log'
    
    lun = -1
    
    if keyword_set(test) then begin
        datdir = locdatdir
        figdir = locfigdir
        load = 1
        movefile = 0
        combine = 0
        logxls = 0
    endif
    
    if keyword_set(default) then begin
        datdir = syncdatdir
        figdir = syncfigdir
        load = 1
        movefile = 1
        combine = 0
        logxls = 0
    endif
    
    if keyword_set(fullupdate) then begin
        datdir = syncdatdir
        figdir = syncfigdir
        load = 1
        movefile = 1
        combine = 1
        logxls = 1
    endif

    ; update data and plots.
    if keyword_set(load) then begin
        printf, lun, 'Update data and plots ...'
        foreach tid, ids do cusp_save_plot_and_data, tid, /reload, /save_data, test = test
    endif

    filetype = [ $
        'field_and_poynt_freq_band', $
        'mat_spec_fast', $
        'mat_spec_polar', $
        'efluxes']
    newdir = [ $
        'poynt_bands', $
        'fast_bands', $
        'polar_bands', $
        'efluxes']

    ; move files to type based directory.
    if keyword_set(movefile) then begin
        printf, lun, 'Update typed based directories ...'
        for i = 0, n_elements(newdir)-1 do begin
            if keyword_set(fullupdate) then spawn, 'rm -r "'+odir+'/'+newdir[i]+'"'
            cusp_move_files, ids, indir = idir, outdir = odir, $
                newdir = newdir[i], file = filetype[i]
        endfor
    endif
    
    ; combine type based pdfs into one pdf.
    if keyword_set(combine) then begin
        printf, lun, 'Combine pdfs into single pdf ...'
        for i = 0, n_elements(newdir)-1 do begin
            printf, lun, newdir[i]
            tfns = file_search(odir+'/'+newdir[i]+'/*', count = ntfn)
            cmd = 'cd "'+odir+'/'+newdir[i]+'"'
            cmd+= '; /usr/local/bin/gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sPDFSETTING=/prepress -sOutputFile='
            cmd+= '../'+newdir[i]+'.pdf'
            for j = 0, ntfn-1 do cmd+= ' '+file_basename(tfns[j])
            spawn, cmd
        endfor
    endif

    ; update excel form.
    if keyword_set(logxls) then begin
        printf, lun, 'Update log for excel form ...'
        cusp_gen_excel_form, ids, filename = xlsfn, /load, test = test
    endif
    
end



;ids = '1998_0922_23'
;cusp_update_system, ids, /default

;ids = cusp_id('default')
;cusp_update_system, ids, /movefile

ids = cusp_id('all')
cusp_update_system, ids, /fullupdate

end