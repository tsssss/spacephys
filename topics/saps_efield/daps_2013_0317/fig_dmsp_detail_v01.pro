;+
;-

function fig_dmsp_detail_v01, test=test

    id = '2013_0317'
    version = 'v01'
    plot_dir = join_path([googledir(),'works','2024_daps',id,'plot'])
    time_range = time_double(['2013-03-17/09:20','2013-03-17/09:45'])
    ssusi_time = time_double('2013-03-17/09:32')
    ssusi_id = 'energy'
    probe = 'f18'
    prefix = 'dmsp'+probe+'_'

    r1_times = list()
    r2_times = list()

    r1r2_boundary = ['2013-03-17/09:25:50','2013-03-17/09:39:40']
    r2_times.add, ['2013-03-17/09:25:06', r1r2_boundary[0]]
    r1_times.add, [r1r2_boundary[0], '2013-03-17/09:27:32']
    
    r1_times.add, ['2013-03-17/09:38:40', r1r2_boundary[1]]
    r2_times.add, [r1r2_boundary[1], '2013-03-17/09:42:52']
    
    r1_color = sgcolor('wheat')    
    r2_color = sgcolor('lavender')

    db_xyz_var = dmsp_read_bfield_madrigal(time_range, probe=probe)
    b0_xyz_var = dmsp_read_bfield_madrigal(time_range, probe=probe, read_b0=1)
    v_var = dmsp_read_ion_vel_madrigal(time_range, probe=probe)
    dmsp_mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id=ssusi_id)
    mlat_vars = dmsp_read_mlat_vars(time_range, probe=probe, errmsg=errmsg)
    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]

    ele_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='e', errmsg=errmsg)
    ion_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='p', errmsg=errmsg)
    ele_eflux_var = dmsp_read_eflux(time_range, probe=probe, species='e', errmsg=errmsg)
    
    b_vec = get_var_data(db_xyz_var, times=times)
    index = where(finite(snorm(b_vec)), count)
    b_vec = sinterpol(b_vec[index,*], times[index], times)
    store_data, db_xyz_var, times, b_vec
    
    b_vec = get_var_data(b0_xyz_var, times=times)
    bmag_var = prefix+'bmag'
    store_data, bmag_var, times, snorm(b_vec)
    add_setting, bmag_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', '|B|', $
        'unit', 'nT')
    options, bmag_var, 'ytitle', '|B| (nT)'

    ; convert dB from xyz to fac.
    fac_labels = [tex2str('perp')+','+['out','west'],tex2str('parallel')]
    db_xyz = get_var_data(db_xyz_var, times=times)
    db_fac = db_xyz
    fmlt = get_var_data(mlt_var, at=times)
    fmlat = get_var_data(mlat_var, at=times)
    theta = (fmlt*15-90)*constant('rad')
    r_hat = [[cos(theta)],[sin(theta)]]
    w_hat = [[cos(theta+0.5*!dpi)],[sin(theta+0.5*!dpi)]]
    db_fac[*,0] = r_hat[*,0]*db_xyz[*,0]+r_hat[*,1]*db_xyz[*,1]
    db_fac[*,1] = w_hat[*,0]*db_xyz[*,0]+w_hat[*,1]*db_xyz[*,1]
    db_var = prefix+'db_fac'
    store_data, db_var, times, db_fac
    add_setting, db_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dB', $
        'unit', 'nT', $
        'coord', '', $
        'coord_labels', fac_labels )
;    var = db_var
;    yrange = [-1,1]*500
;    ytickv = [-1,0,1]*300
;    yticks = n_elements(ytickv)-1
;    yminor = 3
;    options, var, 'yrange', yrange
;    options, var, 'ytickv', ytickv
;    options, var, 'yticks', yticks
;    options, var, 'yminor', yminor

    var = prefix+'v_dmsp_xyz'
    options, var, 'labels', 'V'+constant('xyz')
    yrange = [-1.5,5.5]
    ytickv = [-1,2,5]
    yminor = 3
    set_ytick, var, yrange=yrange, ytickv=ytickv, yminor=yminor
    options, var, 'constant', ytickv

    var = prefix+'db_dmsp_xyz'
    options, var, 'labels', 'dB'+constant('xyz')
    yrange = [-1,1]*700
    ytickv = [-1,0,1]*500
    yminor = 5
    set_ytick, var, yrange=yrange, ytickv=ytickv, yminor=yminor
    options, var, 'constant', ytickv

    base = 'fig_dmsp_detail_'+id+'_'+version+'.pdf'
    plot_file = join_path([plot_dir,base])
    print, plot_file
    if keyword_set(test) then plot_file = 0
    fig_size = [6,4]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    uniform_ticklen = -ychsz*0.3*fig_size[1]
    
    margins = [10,6,8,1]
    plot_vars = prefix+['e_en_spec','p_en_spec','v_dmsp_xyz','db_dmsp_xyz']
    nplot_var = n_elements(plot_vars)
    poss = sgcalcpos(nplot_var, margins=margins)
    tplot_options, 'tickinterval', 5*60
    var_labels = prefix+['mlt_plot','mlat','bmag']
    copy_data, prefix+'mlt', prefix+'mlt_plot'
    var = prefix+'mlt_plot'
    mlt = get_var_data(var, times=times)
    index = where(mlt le 0, count)
    if count ne 0 then begin
        mlt[index] += 24d
        store_data, var, times, mlt
    endif
    options, var, 'ytitle', 'MLT (h)'
    var = prefix+'mlat'
    options, var, 'ytitle', 'MLat (deg)'
    
    vars = prefix+['e','p']+'_en_spec'
    options, vars, 'ytitle', 'Energy!C(eV)'

    tpos = poss[*,nplot_var-2]
    tpos[1] = poss[1,nplot_var-1]
    yrange = [0,1]
    set_axis, position=tpos, xrange=time_range, yrange=yrange
    foreach tr, r1_times do begin
        tr = time_double(tr)
        polyfill, tr[[0,1,1,0,0]], yrange[[0,0,1,1,0]], color=r1_color
        tx = mean(tr)
        ty = yrange[1]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*1
        msg = 'R1!C'
        mlt = get_var_data(prefix+'mlt', at=tr[0])
        if mlt le 0 then msg += 'U' else msg += 'D'
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach
    foreach tr, r2_times do begin
        tr = time_double(tr)
        polyfill, tr[[0,1,1,0,0]], yrange[[0,0,1,1,0]], color=r2_color
        tx = mean(tr)
        ty = yrange[1]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*1
        msg = 'R2!C'
        mlt = get_var_data(prefix+'mlt', at=tr[0])
        if mlt ge 0 then msg += 'U' else msg += 'D'
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach
    
    for pid=0,nplot_var-1 do begin
        tpos = poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
    endfor

    tplot, plot_vars, trange=time_range, position=poss, var_label=var_labels, vlab_margin=9, noerase=1
    fig_labels = letters(nplot_var)+') '+['e-','H+','Ion Vel','dB']
    for pid=0,nplot_var-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*8.5
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]
        xyouts, tx,ty,msg, normal=1
        
        if pid eq 0 then begin
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = strupcase('dmsp '+probe)
            xyouts, tx,ty,msg, normal=1, color=sgcolor('white')
        endif
    endfor
    
    tpos = poss[*,-1]
    yrange = [0,1]
    set_axis, position=tpos, xrange=time_range, yrange=yrange
    msg_times = list()
    msg_times.add, minmax(time_double([r2_times[0],r1_times[0]]))
    msg_times.add, minmax(time_double([r2_times[1],r1_times[1]]))
    msgs = ['Dawn','Dusk']
    foreach msg, msgs, tid do begin
        tx = mean(msg_times[tid])
        ty = yrange[1]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*1
        if tid eq 1 then ty = tpos[1]+ychsz*0.5
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach

    if keyword_set(test) then stop
    sgclose
    return, plot_file
    
end

print, fig_dmsp_detail_v01(test=0)
end
