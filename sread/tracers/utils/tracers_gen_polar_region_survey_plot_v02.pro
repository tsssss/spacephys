;+
; v01 is mostly L1
; v02 is for L2.
; 
; Adopted from dmsp_gen_polar_region_survey_plot_v04.
;-

function tracers_gen_polar_region_survey_plot_v02, input_time_range, probe=probe, plot_dir=plot_dir, errmsg=errmsg, test=test, local_root=local_root, tickinterval=tickinterval


    compile_opt idl2
    errmsg = ''
    retval = !null
    file = get_filename()
    version = get_file_version(file)
    

    time_range = time_double(input_time_range)
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'survey_plot','tracers_'+version])

    ; Load data.
    prefix = 'ts'+probe+'_'

    full_time_range = minmax(time_range)+[-1,1]*60
    full_time_range = minmax(time_range)
    mlat_vars = tracers_read_mlat_vars(full_time_range, probe=probe, errmsg=errmsg, update=1)
    if errmsg ne '' then return, retval
    ele_spec_var = tracers_read_en_spec(full_time_range, probe=probe, species='e', errmsg=errmsg)
    if errmsg ne '' then return, retval
    ion_spec_var = tracers_read_en_spec(full_time_range, probe=probe, species='i', errmsg=errmsg)
    if errmsg ne '' then return, retval

    mlat_var = mlat_vars['mlat']
    mlt_var = mlat_vars['mlt']

    ; Load dB.
    db_var = tracers_read_bfield(full_time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval
    ;add_setting, db_var, id='bfield', dictionary('coord','FAC')
    options, db_var, yrange=[-1,1]*600, yminor=5, yticks=2, constant=0

    ; Load E spec.
    de_var = tracers_read_efield(full_time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval


    ; Generate plot.
    if n_elements(plot_dir) eq 0 then plot_dir = join_path([local_root,'%Y','%m%d'])

    plot_files = list()
    the_time_range = time_range
    time = mean(the_time_range)
    mlats = var_get_data(mlat_var, in=the_time_range)
    hem = median(mlats) lt 0 ? 'south' : 'north'
    is_south = median(mlats) lt 0
    hem_str = is_south? 'South': 'North'
        
    margins = [12,6,2,1]
    all_poss = panel_pos(pansize=[1,1]*4, panid=[1,0], xpans=[2,1], ypans=[1], xpad=10, margins=margins, fig_size=fig_size)
    path = apply_time_to_pattern(plot_dir,time)
    base = 'tracers_polar_region_survey_'+strlowcase(hem)+'_'+strjoin(time_string(the_time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_'+probe+'_'+version+'.pdf'
    plot_file = join_path([path,base])
    plot_files.add, plot_file
    if keyword_set(test) then begin
        plot_file = 0
    endif else begin
;        if file_test(plot_file) eq 1 then begin
;            print, plot_file+' exists, skip ...'
;            return, plot_files.toarray()
;        endif
    endelse

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    plot_vars = [ele_spec_var,ion_spec_var,db_var,de_var]
    nplot_var = n_elements(plot_vars)
    plot_labels = letters(nplot_var)+'. '+['e-','Ion','dB','Ex']

    label_vars = [mlat_var,mlt_var]
    options, mlat_var, ytitle='MLat (deg)'
    options, mlt_var, ytitle='MLT (h)'
    vlab_margins = 10

    spec_vars = [ele_spec_var,ion_spec_var]
    options, spec_vars, ytitle='Energy!C(eV)', zticklen=-0.5

    plot_pos = combine_pos(all_poss[*,0,*])
    left_poss = sgcalcpos(nplot_var,position=plot_pos)
    if n_elements(tickinterval) eq 0 then tickinterval = 60
    tplot_options, 'tickinterval', tickinterval
    set_axis, position=plot_pos, xrange=the_time_range, yrange=[0,1]
    tplot, plot_vars, var_label=label_vars, $
        trange=the_time_range, noerase=1, position=left_poss, $
        vlab_margin=vlab_margins
        
    for ii=0,nplot_var-1 do begin
        my_pos = left_poss[*,ii]
        tx = my_pos[0]-xchsz*vlab_margins
        ty = my_pos[3]-ychsz*0.8
        msg = plot_labels[ii]
        xyouts, tx,ty,msg, normal=1
    endfor


;---MLT image.
    big_pos = reform(all_poss[*,1,0])
    
;---Add dB.
    down_sample_cadence = 6
    plot_scale = 10d
    db_color = sgcolor('peru')
    db_unit = 'nT'
    db_scale = 300d     ; nT.
    label_color = sgcolor('black')


    the_color = db_color
    the_unit = db_unit
    the_scale = db_scale
    my_pos = polar_xy_set_axis(big_pos)
    polar_xy_draw_vector, db_var, time_range=the_time_range, $
        mlat_var=mlat_var, mlt_var=mlt_var, $
        data_scale=the_scale, $
        plot_scale=plot_scale, $
        color=the_color, down_sample_cadence=down_sample_cadence
    polar_xy_draw_axis, south=is_south, mlt_label_mlat=55, $
        mlat_tickformat='(A1)'

    ; Add label and scale.
    len = polar_xy_mlat_to_dis(90-plot_scale)
    label_pos = [my_pos[0]+xchsz*3,my_pos[3]-ychsz*1.4]
    tmp = convert_coord(label_pos, normal=1, to_data=1)
    txs = tmp[0]+[0,len]
    tys = tmp[1]+[0,0]
    plots, txs, tys, data=1, color=the_color
    foreach tx, txs do begin
        tmp = convert_coord(tx,tys[0], data=1, to_normal=1)
        ttxs = tmp[0]+[0,0]
        ttys = tmp[1]+[-1,1]*ychsz*0.15
        plots, ttxs, ttys, normal=1, color=the_color
    endforeach
    mid_pos = [mean(txs),tys[0]]
    tmp = convert_coord(mid_pos, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]+ychsz*0.3
    msg = string(the_scale,format='(I0)')+' '+the_unit
    xyouts, tx,ty,msg,normal=1, alignment=0.5, color=the_color


;---Draw orbit.
    line_color = sgcolor('silver')
    mlts = get_var_data(mlt_var, in=the_time_range, times=the_times)
    mlats = get_var_data(mlat_var, in=the_time_range)
    sc_xys = polar_xy_from_mlat_mlt(mlats, mlts)
    sc_xs = sc_xys[*,0]
    sc_ys = sc_xys[*,1]
    minor_times = make_bins(the_time_range, 60, inner=1)
    minor_xys = sinterpol(sc_xys, the_times, minor_times)
    major_times = make_bins(the_time_range, 120, inner=1)
    major_xys = sinterpol(sc_xys, the_times, major_times)
    major_tickns = time_string(major_times,tformat='hh:mm')
    major_xxs = major_xys[*,0]
    major_yys = major_xys[*,1]

    my_pos = big_pos
    my_pos = polar_xy_set_axis(my_pos)

    oplot, sc_xs, sc_ys, color=line_color

    set_circ, fill=1
    plots, minor_xys[*,0], minor_xys[*,1], psym=8, symsize=0.5, color=line_color

    foreach msg, major_tickns, ii do begin
        tmp = convert_coord(major_xxs[ii],major_yys[ii], data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]+ychsz*0.4
        plots, tmp[0], tmp[1], normal=1, psym=8, symsize=0.5, color=label_color
        if ii eq 1 then tx = tmp[0]-xchsz*3
        xyouts, tx,ty,normal=1, msg, alignment=0.5, color=label_color
    endforeach

    ; Add letters.
    plot_labels = letters([0,1]+nplot_var)+'. '
    my_pos = big_pos
    tx = my_pos[0]+xchsz*0.5
    ty = my_pos[3]-ychsz*1.0
    msg = plot_labels
    xyouts, tx,ty,msg, normal=1

    if keyword_set(test) then stop
    sgclose

    
    return, plot_files.toarray()

end


compile_opt idl2
test = 0
time_range = ['2026-02-16/04:00','2026-02-16/04:05']
time_range = ['2026-02-16/04:00','2026-02-16/04:18']
probe = '2'

tr = ['2026-02-16','2026-02-17']
mlat_vars = tracers_read_mlat_vars(tr, probe=probe, errmsg=errmsg, update=1)
mlat_var = mlat_vars['mlat']
mlats = var_get_data(mlat_var, times=times)
min_mlat = 50d
index = where(mlats ge min_mlat, count)
north_trs = times[time_to_range(index)]
nnorth_tr = n_elements(north_trs[*,0])
index = where(mlats le -min_mlat, count)
south_trs = times[time_to_range(index)]
nsouth_tr = n_elements(south_trs[*,0])

; Load data first.
full_time_range = minmax(tr)
mlat_vars = tracers_read_mlat_vars(full_time_range, probe=probe, errmsg=errmsg, update=1)
ele_spec_var = tracers_read_en_spec(full_time_range, probe=probe, species='e', errmsg=errmsg)
ion_spec_var = tracers_read_en_spec(full_time_range, probe=probe, species='i', errmsg=errmsg)
db_var = tracers_read_bfield(full_time_range, probe=probe, errmsg=errmsg)
de_var = tracers_read_efield(full_time_range, probe=probe, errmsg=errmsg)

; for ii=0,nnorth_tr-1 do begin
;     the_time_range = reform(north_trs[ii,*])
;     ; Check if it's a full orbit.
;     the_mlats = var_get_data(mlat_var, at=the_time_range)
;     if abs(mean(abs(the_mlats))-min_mlat) ge 2 then begin
;         print, 'Not a full orbit, skip ...'
;         continue
;     endif
; 
;     files = tracers_gen_polar_region_survey_plot_v02(the_time_range, probe=probe, test=test)
;     print, files
; endfor

for ii=0,nsouth_tr-1 do begin
    the_time_range = reform(south_trs[ii,*])
    ; Check if it's a full orbit.
    the_mlats = var_get_data(mlat_var, at=the_time_range)
    if abs(mean(abs(the_mlats))-min_mlat) ge 2 then begin
        print, 'Not a full orbit, skip ...'
        continue
    endif

    files = tracers_gen_polar_region_survey_plot_v02(the_time_range, probe=probe, test=test)
    print, files
endfor

stop
files = tracers_gen_polar_region_survey_plot_v02(time_range, probe=probe, test=test, tickinterval=60)
print, files

end
