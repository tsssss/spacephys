;+
; An overview figure, showing the conjunction of wave and arc, and ion beams.
;-

function fig_2013_0501_field_v02, plot_file, event_info=event_info

test = 1



;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
    plasma_param = event_info['plasma_param']

;---Settings.
    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']
    snapshot_time = event_info['snapshot_time']
    fac_labels = event_info['fac_labels']
    bar_times = make_bins(time_range,600, inner=1)
    
    
    xticklen_chsz = -0.2
    yticklen_chsz = -0.35
    
    ; n var.
    n_var = prefix+'plot_density'
    get_data, prefix+'density_hope', times, data, limits=lim
    index = where(finite(data))
    data = interpol(data[index], times[index], times)
    store_data, n_var, times, data, limits=lim
    
    var = n_var
    yrange = [0.1,3]
    ytickv = [0.1,1]
    yticks = n_elements(ytickv)-1
    yminor = 10
    constant = [1]
    options, var, 'ystyle', 1
    options, var, 'ylog', 1
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant


    ; vexb var.
    vexb_var = prefix+'vexbdot0_fac'
    vexb_vars = prefix+'vexbdot0_'+['para','west','out']
    vexb_labels = 'ExB V!D'+fac_labels+'!N'
    stplot_split, vexb_var, newnames=vexb_vars, labels=vexb_labels
    
    var = vexb_vars[1]
    yrange = [-1,1]*600
    ystep = 400
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 4
    constant = [0,-1,1]*ystep
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant
    
    var = vexb_vars[2]
    yrange = [-1,1]*300
    ystep = 200
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 2
    constant = [0,-1,1]*ystep
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant
    
    foreach var, vexb_vars do begin
        uniform_time, var, 5
    endforeach


    ; e vars.
    e_var = prefix+'edot0_fac'
    e_vars = prefix+'edot0_'+['para','west','out']
    e_labels = 'E!D'+fac_labels+'!N'
    stplot_split, e_var, newnames=e_vars, labels=e_labels
    
    var = e_vars
    yrange = [-1,1]*190
    ystep = 100
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    constant = [0,-1,1]*ystep
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant
    
    
    ; dB vars.
    db_var = prefix+'b1_fac'
    var = db_var
    yrange = [-1,1]*60
    ystep = 50
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', 0

    ; B var.
    b_var = prefix+'b_sm'
    var = prefix+'b_gsm'
    get_data, var, times, b_gsm, limits=lim
    b_sm = cotran(b_gsm, times, 'gsm2sm')
    store_data, b_var, times, b_sm
    add_setting, b_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', 'nT', $
        'coord', 'SM', $
        'coord_labels', constant('xyz') )
    yrange = [-90,190]
    ystep = 100
    yminor = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    constant = ytickv[lazy_where(ytickv,'()', yrange)]
    var = b_var
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant


    ; O+ outflow.
    ;vars = rbsp_read_en_spec_combo(time_range, probe=probe, species='o')
    spec_vars = prefix+[$
        'o_en_spec','o_en_spec_'+['anti','para'], $
        'p_en_spec','p_en_spec_'+['anti','para'] ]
    foreach var, spec_vars do begin
        get_data, var, times, data, vals
        index = where(data eq 0 or finite(data,nan=1), count)
        if count ne 0 then begin
            data[index] = 1e-2
            store_data, var, times, data, vals
        endif
    endforeach
    yrange = [5,5e4]
    log_yrange = alog10(yrange)
    log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    yminor = 10
    foreach tx, ytickv, ii do begin
        if tx eq 1 then begin
            ytickn[ii] = '1'
        endif else if tx eq 10 then begin
            ytickn[ii] = '10'
        endif
    endforeach

    var = [spec_vars,prefix+['p_en_spec','e_en_spec']]
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', ytickv
    options, var, 'ytitle', 'Energy!C(eV)'
    options, var, 'zticklen', -0.3
    options, var, 'zminor', 10
    
    var = prefix+['p_en_spec','p_en_spec_'+['anti','para']]
    options, var, 'zrange', [1e4,1e6]
    var = prefix+['o_en_spec','o_en_spec_'+['anti','para']]
    options, var, 'zrange', [5e3,5e5]    
    
    var = prefix+'e_en_spec'
    yrange = [15,5e4]
    log_yrange = alog10(yrange)
    log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    yminor = 10
    foreach tx, ytickv, ii do begin
        if tx eq 1 then begin
            ytickn[ii] = '1'
        endif else if tx eq 10 then begin
            ytickn[ii] = '10'
        endif
    endforeach
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', ytickv
    options, var, 'ytitle', 'Energy!C(eV)'
    zrange = [1e5,1e10]
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,2,inner=1)
    ztickv = 10.^log_ztickv
    zticks = n_elements(ztickv)-1
    options, var, 'zrange', zrange
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    options, var, 'zcharsize', label_size
    
    
    var = prefix+['p_en_spec','p_en_spec_'+['anti','perp','para']]
    options, var, 'color_table', 63

    var = prefix+['o_en_spec','o_en_spec_'+['anti','perp','para']]
    options, var, 'color_table', 64
    
    var = prefix+'e_en_spec'
    options, var, 'color_table', 65





    
    
;---Plot.
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_field_v02.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz

    plot_vars = [n_var,e_vars[1:2],db_var,b_var,vexb_vars[1:2]]
    nplot_var = n_elements(plot_vars)
    nvar = nplot_var
    fig_letters = letters(nvar*2)
    fig_labels = fig_letters[0:nplot_var-1]+') '+['N','E!D'+fac_labels[1:2],$
        'dB','SM B','V!D'+fac_labels[1:2]]
        

    margins = [10,4,10,1]
    label_size = 0.8

    tmp = panel_pos(plot_file, pansize=[3,0.75], $
        xpans=[1,1], ypans=1+fltarr(nplot_var), margins=margins, $
        fig_size=fig_size, xpad=17)

    poss_left = reform(tmp[*,0,*])
    poss_right = reform(tmp[*,1,*])

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    

;---tplot panels.
    plot_poss = poss_left

    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, plot_vars[ii], 'xticklen', xticklen
        options, plot_vars[ii], 'yticklen', yticklen
    endfor

    tplot, plot_vars, trange=time_range, position=plot_poss, noerase=1
    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        fig_label = fig_labels[ii]
        tx = 2*xchsz
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,fig_label, normal=1
    endfor
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    timebar, snapshot_time, color=sgcolor('red')

    
    spec_poss = poss_right
    spec_vars = prefix+['e_en_spec',$
        'p_en_spec', 'p_en_spec_'+['para','anti'], $
        'o_en_spec', 'o_en_spec_'+['para','anti'] ]
    spec_labels = fig_letters[nplot_var:-1]+') '+['e-',$
        'H+','H+ '+['para','anti'], $
        'O+','O+ '+['para','anti'] ]
    
    for ii=0,nplot_var-1 do begin
        tpos = spec_poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, spec_vars[ii], 'xticklen', xticklen
        options, spec_vars[ii], 'yticklen', yticklen
    endfor
    
    foreach var, prefix+['p','o']+'_en_spec' do begin
        vid = where(spec_vars eq var)
        tpos = spec_poss[*,vid+2]
        cbpos = spec_poss[*,vid]
        cbpos[0] = tpos[2]+xchsz*1
        cbpos[2] = cbpos[0]+xchsz*1
        cbpos[1] = tpos[1]
        tpos = spec_poss[*,vid]
        cbpos[3] = tpos[3]
        options, var, 'zposition', cbpos
        options, var, 'no_color_scale', 0
        vars = var+'_'+['para','anti']
        options, vars, 'no_color_scale', 1
    endforeach
    
    tplot, spec_vars, trange=time_range, position=spec_poss, noerase=1
    for ii=0,nplot_var-1 do begin
        tpos = spec_poss[*,ii]
        fig_label = spec_labels[ii]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,fig_label, normal=1
    endfor
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    timebar, snapshot_time, color=sgcolor('red')


;---Trace O+ data.
    model_time = time_double('2013-05-01/07:38:03')
    beam_dis = 2.1
    color_beam = sgcolor('deep_pink')
    color_conics = sgcolor(['green','blue'])
    test_pitch_angles = [175d,162]


    test_info_list = list()
    test_info_list.add, dictionary($
        'species', 'o', $
        'trs', time_double('2013-05-01/'+['07:40:53','07:43:09']), $
        'ens', [6200,1300] )
;    test_info_list.add, dictionary($
;        'species', 'p', $
;        'trs', time_double('2013-05-01/'+['07:39:22','07:40:53']), $
;        'ens', [6000,300] )

    trace_input_list = list()
    foreach info, test_info_list do begin
        species = info['species']
        trs = info['trs']
        ens = info['ens']
        
        test_time_range = trs
        
    ;---To get the trace info.
        get_data, prefix+species+'_en_spec', times, data
        index = lazy_where(times, '[]', trs)
        test_times = times[index]
        log_ens = alog10(ens)
        foreach time, test_times, time_id do begin
            energy = 10.^(log_ens[0]+(log_ens[1]-log_ens[0])/(trs[1]-trs[0])*(time-trs[0]))
            
            trace_input_list.add, dictionary($
                'species', species, $
                'time', time, $
                'energy', energy, $
                'pitch_angles', test_pitch_angles, $
                'beam_dis', beam_dis, $
;                'mod_time', time, $
                'mod_time', model_time, $
                'model', 't89', $
                'igrf', 0 )
        endforeach
    endforeach
    
    
    
    ; Collect plot_info for input data.
    plot_info = dictionary()
    foreach info, trace_input_list do begin
        species = info['species']
        
        if ~plot_info.haskey(species) then begin
            plot_info[species] = dictionary($
                'times', list(), $
                'energys', list() )
        endif
        
        (plot_info[species])['times'].add, info['time']
        (plot_info[species])['energys'].add, info['energy']        
    endforeach

    ; Plot plot_info.
    foreach species, plot_info.keys() do begin
        north_var = prefix+species+'_en_spec_anti'
        south_var = prefix+species+'_en_spec_para'
        
        pid = where(spec_vars eq north_var)
        tpos = spec_poss[*,pid]
        xrange = time_range
        yrange = get_setting(north_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, ylog=1, $
            nodata=1, noerase=1, position=tpos

        times = ((plot_info[species])['times']).toarray()
        energys = ((plot_info[species])['energys']).toarray()
        plots, times, energys, data=1, psym=1, symsize=0.2, color=sgcolor('black')
    endforeach

    

;---Trace and plot output.
    r_gsm_var = prefix+'r_gsm'
    trace_output_list = list()
    foreach info, trace_input_list do begin        
        time = time_double(info['time'])
        info['r_gsm'] = get_var_data(r_gsm_var, at=time)
        info['time'] = time
        par_var = info['model']+'_var'
        info['par'] = reform(get_var_data(par_var, at=time))
    
        trace_output = trace_ion_to_ionosphere(time, _extra=info.tostruct())
        trace_output['species'] = info['species']
        trace_output['energy'] = info['energy']
        trace_output_list.add, trace_output
    endforeach
    
    ; Collect plot_info.
    plot_info = dictionary()
    dirs = ['anti','para']
    foreach info, trace_output_list do begin
        species = info['species']
        foreach dir, dirs do begin
            the_var = prefix+species+'_en_spec_'+dir
            if ~plot_info.haskey(the_var) then begin
                plot_info[the_var] = dictionary($
                    'conic_times', list(), $
                    'beam_times', list(), $
                    'energys', list(), $
                    'conic_dis', list() )
            endif
            
            key = (dir eq 'anti')? 'north': 'south'
            the_output = info[key]
            (plot_info[the_var])['conic_times'].add, the_output['conic_time']
            (plot_info[the_var])['beam_times'].add, the_output['beam_time']
            (plot_info[the_var])['energys'].add, info['energy']
            (plot_info[the_var])['conic_dis'].add, the_output['conic_dis']
        endforeach
    endforeach
    
    
    ; Plot plot_info.
    foreach var, plot_info.keys() do begin
        pid = where(spec_vars eq var)
        tpos = spec_poss[*,pid]
        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, ylog=1, $
            nodata=1, noerase=1, position=tpos

        the_output = plot_info[var]
        conic_times = (the_output['conic_times']).toarray()
        beam_times = (the_output['beam_times']).toarray()
        energys = (the_output['energys']).toarray()

        foreach color_conic, color_conics, conic_id do begin        
            plots, conic_times[*,conic_id], energys, data=1, color=color_conic, psym=1, symsize=0.2
        endforeach
        
        plots, beam_times, energys, data=1, color=color_beam, psym=1, symsize=0.2
    endforeach
    
    
    var = prefix+'o_en_spec_para'
    pid = where(spec_vars eq var)
    tpos = spec_poss[*,pid]
    xrange = time_range
    yrange = get_setting(var, 'yrange')
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, ylog=1, $
        nodata=1, noerase=1, position=tpos
    tx = snapshot_time
    ty = yrange[0]
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]+xchsz*0.5
    tx = tpos[0]+xchsz*0.5
    the_info = plot_info[var]
    conic_dis = the_info['conic_dis'].toarray()
    foreach color_conic, color_conics, conic_id do begin
        ty = tpos[1]+(0.2+conic_id)*ychsz
        the_dis = mean(conic_dis[*,conic_id])
        msg = 'Conic, PA '+string(test_pitch_angles[conic_id],format='(I0)')+' deg, '+string(the_dis,format='(F3.1)')+' Re'
        xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color_conic
    endforeach
    
    ty = tpos[1]+(0.2+2)*ychsz
    msg = 'Beam, '+string(beam_dis,format='(F3.1)')+' Re'
    xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color_beam
    
    
;    foreach trace_output, trace_output_list do begin
;        foreach key, ['north','south'] do begin
;            the_output = trace_output[key]
;            print, the_output['conic_dis'], the_output['beam_dis']
;        endforeach
;    endforeach


    
    
    if keyword_set(test) then stop
    sgclose
    
    

end

print, fig_2013_0501_field_v02(event_info=event_info)
end