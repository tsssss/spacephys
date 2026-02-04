
function micro_injection_themis_2009_0619_fig_themis_v02, test=test

    plot_dir = srootdir()
    base = 'micro_injection_themis_2009_0619_fig_themis_v02.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0

    margins = [14,4,8,1]
    poss = panel_pos(plot_file, fig_size=fig_size, margins=margins)
    fig_size = [6.5,5]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(1,margins=margins)

    probe = 'b'
    prefix = 'th'+probe+'_'
    time_range = time_double('2009-06-19/'+['00:00','05:30'])

    ; Load data.
    mission_probe = 'th'+probe

    ; B field and ion velocity.
    b_gsm_var = themis_read_bfield(time_range, probe=probe, errmsg=errmsg, id='fgs')
    u_gsm_var = themis_read_ion_vel(time_range, probe=probe, errmsg=errmsg, id='peir')    

    ; SST.
    datatype = 'psef'
    prefix2 = prefix+datatype+'_'
    en_high_var = prefix2+'en_eflux'
    pa_high_var = prefix2+'an_eflux_pa'
    if check_if_update(en_high_var, time_range) then begin
        thm_part_load, data_type=datatype, probe=probe, trange=time_range
        thm_part_getspec, data_type=datatype, probe=probe, trange=time_range, outputs='energy'
        thm_part_getspec, data_type=datatype, probe=probe, trange=time_range, outputs='pa'
        options, en_high_var, requested_time_range=time_range
        options, pa_high_var, requested_time_range=time_range
    endif

    unit = 'eV/cm!E2!N-s-sr-eV'
    zrange = [1e2,1e7]
    ztickv = [1e2,1e3,1e4,1e5,1e6,1e7]
    zrange = [1e1,1e5]
    ztickv = [1e1,1e2,1e3,1e4,1e5]
    ztickv_log = alog10(ztickv)
    zticks = n_elements(ztickv)-1
    ztickn = '10!U'+string(ztickv_log,format='(I0)')
    ztickn[0:*:2] = ' '
    options, [en_high_var], color_table=40, no_interp=1, $
        ytitle='Energy!C(eV)', ztitle=unit, $
        zrange=zrange, zstyle=1, zlog=1, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=9, $
        yrange=[3.1e4,7.2e5], ystyle=1, ylog=1, ytickv=[5e4,5e5], ytickname='10!U'+['4','5'], yticks=1, yminor=9
    options, [pa_high_var], color_table=40, no_interp=1, $
        ytitle='PA!C(deg)', ztitle=unit, $
        zrange=zrange, zstyle=1, zlog=1, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=9, $
        yrange=[0,180], ystyle=1, ylog=0, ytickv=[30,90,150], ytickname=['30','90','150'], yticks=2, yminor=6
    foreach var, [en_high_var,pa_high_var] do begin
        add_setting, var, smart=1, dictionary($
            'display_type', 'spec' )
    endforeach
    
    ; ESA.
    datatype = 'peef'
    prefix2 = prefix+datatype+'_'
    zrange = [1e5,1e8]
    ztickv = [1e5,1e6,1e7,1e8]
    ztickn = '10!U'+['5','6','7','8']
;    ztickn[0:2:*] = ' '
    zticks = n_elements(ztickv)-1
    en_low_var = themis_read_en_spec(time_range, probe=probe, species='e', id='esa_l2')
    options, en_low_var, color_table=40, no_interp=1, $
        zrange=zrange, zstyle=1, zlog=1, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=9, $
        yrange=[1.1e1,2.6e4], ystyle=1, ylog=1, ytickv=[1e2,1e3,1e4], ytickname='10!U'+['2','3','4'], yticks=2, yminor=9
    
    plot_vars = [en_high_var,pa_high_var,en_low_var,b_gsm_var, u_gsm_var]
    nplot_var = n_elements(plot_vars)
    panel_labels = letters([0,nplot_var]+1) +'. '+['e- EN high','e- PA high','e- EN low','B GSM','V!Dion!N GSM']
    panel_labels = letters([0,nplot_var]) +'. '+['e- EN high','e- PA high','e- EN low','B GSM','V!Dion!N GSM']
    right_pos = poss
    panel_poss = sgcalcpos(nplot_var, margins=[0,0,0,0], region=right_pos)

    uniform_ticklen = -ychsz*fig_size[0]*0.15
    for pid=0,nplot_var-1 do begin
        tpos = panel_poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        var = plot_vars[pid]
        options, var, xticklen=xticklen, yticklen=yticklen
        is_spec = var_get_setting(var, 'spec')
        if is_spec then begin
            zticklen = -0.5
            options, var, zticklen=zticklen
        endif
    endfor

    options, b_gsm_var, yrange=[-1,1]*25, ytickv=[-1,0,1]*20, yticks=2, yminor=5, constant=[0]
    options, u_gsm_var, yrange=[-1,1]*220, ytickv=[-1,0,1]*200, yticks=2, yminor=5, constant=[0]

    tplot_options, 'tickinterval', 3600d
    tplot, plot_vars, trange=time_range, position=panel_poss, noerase=1
    for pid=0,nplot_var-1 do begin
        tpos = panel_poss[*,pid]
        tx = tpos[0]-xchsz*12
        ty = tpos[3]-ychsz*0.8
        msg = panel_labels[pid]
        xyouts, tx,ty,msg, normal=1
    endfor

    ; Add label.
    tpos = panel_poss[*,0]
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'TH-'+strupcase(probe)/"++++++++++++++++++
    xyouts, tx,ty,msg, normal=1, color=sgcolor('white')

    ; Add microinjection times.
    mi_times = time_double([ $
        '2009-06-19/03:27:59', $
        '2009-06-19/03:14:00', $
        '2009-06-19/03:45:30' ])
    timebar, mi_times, color=sgcolor('red'), linestyle=2

    if keyword_set(test) then stop
    sgclose

    return, plot_file

end


test = 0
print, micro_injection_themis_2009_0619_fig_themis_v02(test=test)
end
