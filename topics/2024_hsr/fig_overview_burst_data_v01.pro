
function fig_overview_burst_data_v01, test=test

    version = 'v01'
    time_range = time_double(['2013-06-07/00:00','2013-06-07/09:00'])
    probes = ['a','b']
    plot_time_range = time_double(['2013-06-07/00:30','2013-06-07/08:30'])
    fc_time = '2013-06-07/01:10'
    fc_time = !null
    

;---Load AE and Dst.
    dst_var = 'dst'
    if check_if_update(dst_var, time_range) then begin
        omni_read_index, time_range
        ystep = 50
        yrange = minmax([0,minmax(get_var_data(dst_var))])
        yrange += [-1,1]*abs(total(yrange*[-1,1]))*0.05
        set_ytick, dst_var, ystep=ystep, yrange=yrange
        options, dst_var, 'constant', [-50,0]
        ae_var = 'ae'
        yrange = [0,max(get_var_data(ae_var))*1.05]
        ystep = 600
        set_ytick, ae_var, ystep=ystep, yrange=yrange
    endif



;---Load RBSP data.
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'

        r_gsm_var = rbsp_read_orbit(time_range, probe=probe)
        b_gsm_var = rbsp_read_bfield(time_range, probe=probe)
        e_spec_var = rbsp_read_efield_spec(time_range, probe=probe)
        mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)


    ;---Model B field.
        external_models = ['t89','t96','t01','t04s']
        nmodel = n_elements(external_models)
        bmod_vars = list()
        foreach external_model, external_models do begin
            bmod_vars.add, lets_read_geopack_bfield(orbit_var=r_gsm_var,internal_model='igrf',external_model=external_model)
        endforeach
        bmod_vars = bmod_vars.toarray()
        bmod_var = bmod_vars[0]
        db_gsm_vars = lets_decompose_bfield(b_var=b_gsm_var, bmod_var=bmod_var)

    ;---B tilt.
        perigee_dis = 2.
        b_tilt_vars = list()
        foreach var, [b_gsm_var,bmod_vars] do begin
            var_out = streplace(var,'gsm','tilt')
            b_tilt_vars.add, lets_calc_vec_elev(var, coord='sm', var_info=var_out)
        endforeach
        b_tilt_vars = b_tilt_vars.toarray()
        
        b_tilt_var = b_tilt_vars[0]
        bmod_tilt_vars = b_tilt_vars[1:*]
        db_tilt_vars = bmod_tilt_vars
        foreach var_in, db_tilt_vars, vid do begin
            var_out = streplace(var_in, 'b_tilt','db_tilt')
            db_tilt_vars[vid] = lets_subtract_vars(b_tilt_var,var_in,save_to=var_out)
            yrange = [min(get_var_data(var_out))*1.05,20]
            set_ytick, var_out, ystep=25, yrange=yrange
            options, var_out, 'constant', 0
        endforeach
        
        labels = ['Obs',strupcase(external_models)]
        colors = sgcolor(['red','blue','green','orange','purple'])
        b_tilt_combo_var = prefix+'b_tilt_combo'
        b_tilt_combo_var = stplot_merge(b_tilt_vars, output=b_tilt_combo_var, labels=labels[0:nmodel], colors=colors[0:nmodel])
        set_ytick, b_tilt_combo_var, ystep=30, yrange=[0,90]
        options, b_tilt_combo_var, ytitle='(nT)'
        
        db_tilt_combo_var = prefix+'db_tilt_combo'
        db_tilt_combo_var = stplot_merge(db_tilt_vars, output=db_tilt_combo_var, labels=labels[1:nmodel], colors=colors[1:nmodel])
        set_ytick, db_tilt_combo_var, ystep=30, yrange=[0,90]
        options, db_tilt_combo_var, ytitle='(nT)'


    ;---B mag.
        b_mag_vars = list()
        foreach var, [b_gsm_var,bmod_vars] do begin
            var_out = streplace(var,'gsm','mag')
            b_mag_vars.add, lets_calc_vec_mag(var, var_info=var_out)
        endforeach
        b_mag_vars = b_mag_vars.toarray()
    
        b_mag_var = b_mag_vars[0]
        bmod_mag_vars = b_mag_vars[1:*]
        db_mag_vars = bmod_mag_vars
        foreach var_in, db_mag_vars, vid do begin
            var_out = streplace(var_in, 'b_mag','db_mag')
            db_mag_vars[vid] = lets_subtract_vars(b_mag_var,var_in,save_to=var_out)
        endforeach
    
        labels = ['Obs',strupcase(external_models)]
        colors = sgcolor(['red','blue','green','orange','purple'])
        ;b_mag_combo_var = prefix+'b_mag_combo'
        ;b_mag_combo_var = stplot_merge(b_mag_vars, output=b_mag_combo_var, labels=labels, colors=colors)
    
        db_mag_combo_var = prefix+'db_mag_combo'
        db_mag_combo_var = stplot_merge(db_mag_vars, output=db_mag_combo_var, labels=labels[1:nmodel], colors=colors[1:nmodel])
    

    ;---Spec combo.
        fc_vars = list()
        foreach species, ['e','p'] do fc_vars.add, rbsp_read_gyro_freq(time_range, probe=probe, species=species)
        var = prefix+'fce_half'
        fce = get_var_data(prefix+'fce', times=times)
        options, prefix+'fce', labels='f!Dc,e'
        options, prefix+'fcp', labels='f!Dc,H'
        ;options, prefix+'fco', labels='f!Dc,O'
        ;options, prefix+'fche', labels='f!Dc,He'

        ;var = prefix+'flh'
        ;fcp = get_var_data(prefix+'fcp', times=times)
        ;store_data, var, times, fcp*43, limits={labels:'f!DLH!N'}
        ;fc_vars.add, var
        fc_vars = fc_vars.toarray()
        fc_colors = get_color(n_elements(fc_vars))
        foreach var, fc_vars, ii do options, var, 'colors', fc_colors[ii]

        e_spec_combo_var = e_spec_var+'_combo'
        store_data, e_spec_combo_var, data=[e_spec_var,fc_vars]
        options, e_spec_combo_var, 'yrange', get_setting(e_spec_var,'yrange')
        options, e_spec_combo_var, 'labels', ''
    endforeach


;---Assemble plot_vars.
    plot_vars = list()
    fig_labels = list()
    ;plot_vars.add, 'dst'
    ;plot_vars.add, 'ae'

    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        plot_vars.add, prefix+'e_spec_combo'
        fig_labels.add, strupcase('rbsp-'+probe+'')
    endforeach

    ; Some combined vars.
    labels = strupcase('rbsp-'+probes)
    colors = sgcolor(['red','blue'])

    var = 'b_tilt'
    vars = 'rbsp'+probes+'_b_tilt'
    var = stplot_merge(vars, output=var, colors=colors, labels=labels)
    yrange = [0,40]
    set_ytick, var, ystep=20, yrange=yrange, yminor=2
    options, var, ytitle='(deg)', constant=[10,20,30]
    plot_vars.add, var
    fig_labels.add, 'B elev'
    

    var = 'dis'
    vars = 'rbsp'+probes+'_dis'
    var = stplot_merge(vars, output=var, colors=colors, labels=labels)
    yrange = [1,6]
    set_ytick, var, ytickv=[1,3,5], yrange=yrange, yminor=2
    options, var, ytitle='(Re)', constant=[1,3,5]
    plot_vars.add, var
    fig_labels.add, '|R|'
    

    var = 'mlt'
    vars = 'rbsp'+probes+'_mlt'
    var = stplot_merge(vars, output=var, colors=colors, labels=labels)
    yrange = [-1,1]*10
    set_ytick, var, ystep=6, yrange=yrange, yminor=3
    options, var, ytitle='(Re)', constant=[-1,0,1]*6
    plot_vars.add, var
    fig_labels.add, 'MLT'

    plot_vars = plot_vars.toarray()
    nplot_var = n_elements(plot_vars)
    fig_labels = letters(nplot_var)+') '+fig_labels.toarray()
    ypans = [1.5,1.5,0.8,0.5,0.5]
    ypad = [fltarr(nplot_var-1)+0.4]
    
    ; add snapshots for burst data.
    b1_info = orderedhash()
    b1_info['TDS'] = dictionary($
        'time_range', time_double('2013-06-07/04:54'+[':38.000',':38.050']), $
        'probe', 'a' )
    b1_info['Chorus'] = dictionary($
        'time_range', time_double('2013-06-07/06:19'+[':38.000',':38.050']), $
        'probe', 'a' )
    ;    b1_info['ECH'] = dictionary($
    ;        'time_range', time_double('2013-06-07/04:15'+[':20.100',':20.500']), $
    ;        'probe', 'b' )   ; ech
    nburst = n_elements(b1_info.keys())
    ypans = [ypans,1]
    ypad = [ypad,4]
    nvar = nplot_var+1


;---Plot.
    plot_file = join_path([srootdir(),'fig_overview_burst_data_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    if keyword_set(test) then thick = 4 else thick = 20
    fig_size = [6,4]

    pansize = [5,ypans[0]]
    margins = [12,4,8,2]
    all_poss = panel_pos(plot_file, margins=margins, nypan=nvar, $
        ypans=ypans,ypad=ypad, pansize=pansize, fig_size=fig_size)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    poss = all_poss[*,0:nplot_var-1]
    tplot_options, 'tickinterval', 60*60
    tplot_options, 'version', 3
    tplot, plot_vars, trange=plot_time_range, position=poss
    for pid=0,nplot_var-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*(margins[0]-2)
        ty = tpos[3]-ychsz*0.8
        msg = fig_labels[pid]
        xyouts, tx,ty,msg, normal=1
    endfor
    
    ; Add label for fc_vars and burst data.
    nfreq_var = n_elements(fc_vars)
    colors = get_color(nfreq_var)
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        e_spec_var = prefix+'e_spec_combo'
        pid = where(plot_vars eq e_spec_var, count)
        
        if count ne 0 then begin
            tpos = poss[*,pid]
            set_axis, e_spec_var, xrange=plot_time_range, position=tpos, ylog=1
            yrange = get_var_setting(e_spec_var, 'yrange')
            foreach var, fc_vars, var_id do begin
                color = colors[var_id]
                yys = get_var_data(var, times=xxs, settings=settings, in=plot_time_range)
                index = where(finite(yys))
                if n_elements(fc_time) eq 0 then tx = xxs[index[0]] else tx = time_double(fc_time)
                ty = interpol(yys[index], xxs[index], tx)
                ;if product(ty-yrange) ge 0 then continue

                msg = settings.labels
                color = settings.colors

                tmp = convert_coord(tx,ty, data=1, to_normal=1)
                tx = tpos[0]-xchsz*2.5
                ty = (tmp[1]<tpos[3])-ychsz*0.5
                xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color
            endforeach
            
            rbsp_efw_phasef_read_b1_time_rate, plot_time_range, probe=probe
            var = prefix+'efw_vb1_time_rate'
            trs = get_var_data(var, times=uts)
            ntr = n_elements(trs[*,0])
            ;the_trs = time_to_range(uts, pad_times=600)
            color = sgcolor('red')
            for ii=0,ntr-1 do begin
                tr = trs[ii,*]
                txs = tr
                foreach tx,txs,tid do begin
                    tmp = convert_coord(tx,yrange[1], data=1,to_normal=1)
                    txs[tid] = tmp[0]
                endforeach
                tys = tpos[3]-ychsz*0.25
                if max(txs) ge tpos[2] then txs = txs<tpos[2]
                plots, txs,tys, thick=thick, color=color, normal=1
            endfor
        endif
    endforeach
    
    
    ; Add labels for the published event
    var = 'b_tilt'
    pid = where(plot_vars eq var, count)
    event_tr = time_double('2013-06-07/'+['04:52','05:02'])

    if count ne 0 then begin
        tpos = poss[*,pid]
        set_axis, var, xrange=plot_time_range, position=tpos
        yrange = get_var_setting(var, 'yrange')
        
        color = sgcolor('red')
        txs = event_tr
        foreach tx,txs,tid do begin
            tmp = convert_coord(tx,yrange[1], data=1,to_normal=1)
            txs[tid] = tmp[0]
        endforeach
        tys = tpos[3]-ychsz*0.25
        if max(txs) ge tpos[2] then txs = txs<tpos[2]
        plots, txs,tys, thick=thick, color=color, normal=1
    endif
    
    
    b1_vars = list()
    foreach wave_id, b1_info.keys() do begin
        the_info = b1_info[wave_id]
        probe = the_info.probe
        b1_tr = the_info.time_range
        prefix = 'rbsp'+probe+'_'
        be_mgse_var = prefix+'be_mgse_'+wave_id
        data_tr = b1_tr+[-1,1]*60
; del_data, be_mgse_var
        if check_if_update(be_mgse_var, data_tr) then begin
            rbsp_read_burst_efield, data_tr, probe=probe, spin_axis='e'
            be_var = prefix+'be_uvw'
            get_data, be_var, times, be_vec, limits=lim
            ndim = 3
            window = 1d
            width = window/sdatarate(times)
            index = where(abs(be_vec[*,2]) ge 4e2, count)
            ntime = n_elements(times)
            pad_nrec = 0.05/sdatarate(times)
            if count ne 0 then begin
                bad_trs = time_to_range(index, time_step=1)
                nbad_tr = n_elements(bad_trs[*,0])
                for ii=0,nbad_tr-1 do begin
                    i0 = (bad_trs[ii,0]-pad_nrec)>0
                    i1 = (bad_trs[ii,1]+pad_nrec)<(ntime-1)
                    be_vec[i0:i1,2] = !values.f_nan
                endfor
            endif
            for ii=0,ndim-1 do be_vec[*,ii] -= smooth(be_vec[*,ii],width,nan=1,edge_zero=1)
            index = where(finite(be_vec[*,2],nan=1), count)
            if count ne 0 then be_vec[index,2] = 0
            be_mgse = cotran_pro(be_vec, times, coord_msg=['rbsp_uvw','rbsp_mgse'], probe=probe)
            store_data, be_mgse_var, times, be_mgse, limits=lim
            add_setting, be_mgse_var, smart=1, dictionary($
                'display_type', 'vector', $
                'short_name', 'E', $
                'unit', 'mV/m', $
                'coord', 'rbsp_mgse' )
            options, be_mgse_var, labels='MGSE E!D'+constant('xyz')+'!N', requested_time_range=data_tr
        endif
        b1_vars.add, be_mgse_var
    endforeach
    
    
    tplot_options, 'tickinterval', 0.01
    tplot_options, 'version', 1
    b1_vars = b1_vars.toarray()
    nvar = n_elements(b1_vars)
    burst_poss = sgcalcpos(1,nvar, position=all_poss[*,-1], xpad=17)
    burst_letter = (letters(nplot_var+1))[-1]
    foreach wave_id, b1_info.keys(), pid do begin
        tpos = burst_poss[*,pid]
    
        the_info = b1_info[wave_id]
        probe = the_info.probe
        b1_tr = the_info.time_range
        prefix = 'rbsp'+probe+'_'
        be_mgse_var = prefix+'be_mgse_'+wave_id
        data_tr = b1_tr+[-1,1]*60
        tplot, be_mgse_var, trange=b1_tr, position=tpos, noerase=1
        msg = burst_letter+'-'+string(pid+1,format='(I0)')+') '+wave_id
        tx = tpos[0]-xchsz*9
        ty = tpos[3]-ychsz*0.7
        xyouts, tx,ty,msg, normal=1
        
        ; add label to spec.
        if count ne 0 then begin
            ttpos = poss[*,0]
            ttpos[1] = poss[1,-1]
            xrange = plot_time_range
            yrange = [0,1]
            set_axis, position=ttpos, xrange=xrange, yrange=yrange
            foreach tmp, b1_tr, tid do begin
                plots, tmp+[0,0], yrange
                tx = convert_coord(tmp,yrange[0], data=1, to_normal=1)
                tx0 = tx[0]
                ty0 = ttpos[1]
                tx1 = tpos[0]
                ty1 = tpos[3]
                if tid eq 1 then tx1 = tpos[2]
                plots, [tx0,tx1],[ty0,ty1], normal=1
            endforeach
        endif
    endforeach

    if keyword_set(test) then stop
    sgclose
    return, plot_file

end

test = 0
print, fig_overview_burst_data_v01(test=test)
end