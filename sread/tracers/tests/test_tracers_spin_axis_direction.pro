;+
; Test to obtain tracers spin axis direction.
;-

function test_tracers_spin_axis_direction, input_time_range, probe=probe, test=test
    compile_opt idl2

    time_range = time_double(input_time_range)
    tracers_load_spice_kernel, time_range, probe=probe
    plot_dir = srootdir()
    version_str = 'v01'
    base = 'test_tracers_spin_axis_direction_'+time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_'+version_str+'.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0
    fig_size = [12,6]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
  
;---Get the pointing of the axes in TSCS expressed in GSE.
    prefix = 'ts'+probe+'_'
    time_step = 1d ; sec.
    times = make_bins(time_range, time_step)
    ut0 = time_string(times[0],tformat='YYYY-MM-DDThh:mm:ss')
    cspice_str2et, ut0, et0
    ets = et0+times-ut0
    uts = time_string(times,tformat='YYYY-MM-DDThh:mm:ss')
    cspice_str2et, uts, ets
    target = 'TS'+probe
    target_coord = 'tscs'
    sci_id = strupcase(target+'_'+target_coord)
    default_coord = 'gse'
    frame = strupcase(default_coord)
  
    cspice_pxform, sci_id, frame, ets, pxform
    ntime = n_elements(times)
    ndim = 3
    uvw_tscs = dblarr([ndim,ndim])
    uvw_tscs[*,0] = [1,0,0d]
    uvw_tscs[*,1] = [0,1,0d]
    uvw_tscs[*,2] = [0,0,1d]
    uvw_coord = dblarr(ntime,ndim,ndim)
    for ii=0,ntime-1 do begin
        for jj=0,ndim-1 do begin
            uvw_coord[ii,*,jj] = pxform[0:ndim-1,0:ndim-1,ii] ## uvw_tscs[*,jj]
        endfor
    endfor

    plot_vars = list()
    panel_labels = list()
    components = constant('uvw')
    for ii=0,ndim-1 do begin
        var = prefix+components[ii]+'_'+default_coord+'_spice'
        settings = dictionary('coord',default_coord)
        var = var_store(var, uvw_coord[*,*,ii], times, settings=settings)
        if ii ne ndim-1 then begin
            plot_vars.add, var
            panel_labels.add, strupcase(target_coord)+' '+strupcase(components[ii])
            options, var, labels=strupcase(default_coord)+' '+strupcase(components[ii])+constant('xyz')
            options, var, yrange=[-1,1]*1.2, constant=[0]
        endif else begin
            vars = stplot_split(var)
            options, vars, ynozero=1
            sa_vars = vars
            plot_vars.add, sa_vars, extract=1
            panel_labels.add, strupcase(components[ii])+constant('xyz'), extract=1
            labels = strupcase(default_coord)+' '+strupcase(components[ii])+constant('xyz')
            for jj=0,ndim-1 do begin
                options, vars[jj], labels=labels[jj]
            endfor
        endelse
    endfor

;---Remove noise.
    spin_period = 6d  ; sec.
    width = spin_period/time_step
    for jj=0,ndim-1 do begin
        sa_var = sa_vars[jj]
        f0 = var_get_data(sa_var, times=times, settings=settings)
        ; 1. smooth to remove most of the oscillation at spin period.
        f1 = smooth(f0, width, nan=1, edge_mirror=1)
        ; 2. downsample to spin period to fully remove the rest oscilation.
        uts = times[0:*:width]
        nut = n_elements(uts)
        f2 = dblarr(nut-1)
        t2 = dblarr(nut-1)
        for ii=0,nut-2 do begin
            index = where_pro(times,'[]',uts[ii:ii+1])
            f2[ii] = mean(f1[index],nan=1)
            t2[ii] = mean(uts[ii:ii+1])
        endfor
        f3 = interpol(f2, t2, times)

        var = sa_var+'_combo'
        var = var_store(var, [[f0],[f3]], times)
        options, var, colors=sgcolor(['black','red']), labels=['Original','Smoothed']
        index = plot_vars.where(sa_var)
        var = sa_var+'_smooth'
        var = var_store(var, f3, times)
        options, var, labels=settings.labels
        plot_vars[index] = var
    endfor

    r_var = tracers_read_orbit(time_range, probe=probe)
    plot_vars.add, r_var
    panel_labels.add, 'R'

    plot_vars = plot_vars.toarray()
    nplot_var = n_elements(plot_vars)
    panel_labels = letters(nplot_var)+'. '+panel_labels.toarray()



    margins = [12,4,10,2]
    plot_poss = sgcalcpos(nplot_var, margins=margins)

    options, plot_vars, ytitle='(#)'
    tickinterval = 5
    tplot_options, version=3
    tplot_options, 'xtickinterval', tickinterval
    tplot, plot_vars, trange=time_range, vlab_margin=10, position=plot_poss
    for pid=0,nplot_var-1 do begin
        tpos = plot_poss[*,pid]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.7
        msg = panel_labels[pid]
        xyouts, tx,ty,msg, normal=1
    endfor

    tpos = plot_poss[*,0]
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1.0
    msg = 'TS'+strupcase(probe)+' | Predict Kernel'
    xyouts, tx,ty,msg, normal=1

    if keyword_set(test) then stop
    sgclose
    return, plot_file

  return, 1
end

compile_opt idl2
probe = '2'
test = 1
time_range = time_double(['2026-01-01','2026-01-03'])
print, test_tracers_spin_axis_direction(time_range, probe=probe, test=test)
end

