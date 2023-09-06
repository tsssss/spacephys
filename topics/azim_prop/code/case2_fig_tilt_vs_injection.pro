;+
; Plot tilt angle vs electron injection.
;-


;---Settings.
    probes = ['thd','g15','g13','rbspb']

    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data
test = 0

    time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
    project = azim_df_load_project()
    azim_df_load_basic_data, project=project

    energy_range = [50,500]
    foreach probe, probes do begin
        prefix = probe+'_'
        old_var = prefix+'kev_e_flux'
        new_var = prefix+'kev_e_flux_plot'
        get_data, old_var, times, eflux, ebins
        index = where_pro(ebins, energy_range)
        eflux = eflux[*,index]
        ebins = ebins[index]
        nebin = n_elements(ebins)
        if nebin ge 6 then begin
            eflux = eflux[*,0:*:2]
            ebins = ebins[0:*:2]
        endif
        yrange = minmax(eflux)
        ylog_range = alog10(yrange)
        ylog_range = [floor(ylog_range[0]), ceil(ylog_range[1])]
        yrange = 10.^ylog_range
        ylog_tickv = make_bins(ylog_range,1)
        ytickv = 10.^ylog_tickv
        yticks = n_elements(ylog_tickv)-1
        ytickn = '10!U'+string(ylog_tickv,format='(I0)')
        if yticks ge 4 then begin
            ytickn[0:*:2] = ' '
        endif
        yminor = 10
        store_data, new_var, times, eflux, ebins
        add_setting, new_var, /smart, {$
            display_type: 'list', $
            ylog:1, $
            color_table: 52, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'keV', $
            short_name: 'e!U-!N flux', $
            yrange:yrange, $
            yticks:yticks, $
            yminor:yminor, $
            ytickv:ytickv, $
            ytickname:ytickn}
    endforeach


    plot_vars = list()
    foreach probe, probes do begin
        prefix = probe+'_'
        the_var = prefix+'theta'
        options, the_var, 'labels', strupcase(probe)+' tilt'
        options, the_var, 'constant', 0
        plot_vars.add, prefix+['theta','kev_e_flux_plot'], /extract
    endforeach
    plot_vars = plot_vars.toarray()
    nvar = n_elements(plot_vars)

    fig_letters = letters(nvar/2)
    fig_labels = list()
    foreach probe, probes, probe_id do begin
        fig_labels.add, fig_letters[probe_id]+'-'+['1','2']+'.', /extract
    endforeach
    fig_labels = fig_labels.toarray()
    
    
    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'case2_fig_tilt_vs_injection.pdf'])
    
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=5, ysize=8
    margins = [10,4,8,1]
    
    poss = sgcalcpos(nvar, margins=margins, xchsz=xchsz, ychsz=ychsz)
    tplot, plot_vars, trange=time_range, position=poss
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = xchsz*2
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,/normal, fig_labels[ii]
    endfor
    if keyword_set(test) then stop
    sgclose
    


end
