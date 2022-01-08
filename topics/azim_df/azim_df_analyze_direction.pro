;+
; Tell if the given event is along the 4 directions.
;-

function azim_df_analyze_direction_westward, rr
    return, sunitvec(vec_cross(rr,[0.,0,1]))
end
function azim_df_analyze_direction_eastward, rr
    return, -azim_df_analyze_direction_westward(rr)
end
function azim_df_analyze_direction_outward, rr
    return, sunitvec(rr)
end
function azim_df_analyze_direction_earthward, rr
    return, -azim_df_analyze_direction_outward(rr)
end

function azim_df_analyze_direction, candidate, project=project, $
    ; The output.
    angles=the_angles, vmags=the_vmags

    retval = ''
    angle_deviation = 45.
    angle_spread = 20.
    dirs = ['westward','eastward','earthward','outward']
    ndir = n_elements(dirs)

    triad_list = candidate.triad_list
    ntriad = triad_list.length
    angles = fltarr(ntriad,ndir)
    the_vmags = fltarr(ntriad)
    foreach triad, triad_list, triad_id do begin
        the_vmags[triad_id] = triad.vmag_obs_time
        center_r_sm = triad.center_r_sm
        foreach dir, dirs, dir_id do begin
            the_dir = call_function('azim_df_analyze_direction_'+dir, center_r_sm)
            the_angle = atan(the_dir[1],the_dir[0])
            vhat = triad.vhat_obs_time
            angles[triad_id,dir_id] = atan(vhat[1],vhat[0])-the_angle
        endforeach
    endforeach
    angles *= constant('deg')
    index = where(angles le -180, count)
    if count ne 0 then angles[index] += 360
    index = where(angles ge 180, count)
    if count ne 0 then angles[index] -= 360

    lprmsg, ''
    lprmsg, 'Processing candidate '+strjoin(time_string(candidate.time_range),' to ')

    ; Pick out the best direction.
    tmp = min(total(abs(angles),1), index)
    the_angles = reform(angles[*,index])
    the_dir = dirs[index]
    lprmsg, 'The direction: '+the_dir
    lprmsg, 'The angles (deg): '+strjoin(string(the_angles,format='(I0)'),', ')

    ; Check deviation.
    index = where(the_angles gt angle_deviation, count)
    if count ne 0 then begin
        lprmsg, 'Angles do not align to any direction, skip ...'
        return, retval
    endif

    ; Check spread.
    the_angle_spread = stddev(the_angles)
    the_angle_spread = round(the_angle_spread)
    lprmsg, 'Angle spread (deg): '+string(the_angle_spread, format='(I0)')
    if the_angle_spread gt angle_spread then begin
        lprmsg, 'Angle spread is too large, skip ...'
        return, retval
    endif

    vmag_spread = 0.3
    vmag_min = 15.
    ;if candidate.time_range[0]-time_double('2013-10-14/21:55:30') eq 0 then stop
    the_vmag_mean = min(the_vmags)
    the_vmag_mean = round(the_vmag_mean)
    lprmsg, 'vmag mean (km/s): '+string(the_vmag_mean,format='(I0)')
    if the_vmag_mean lt vmag_min then begin
        lprmsg, 'vmag mean is too small, skip ...'
        return, retval
    endif
    the_vmag_spread = stddev(the_vmags)/the_vmag_mean
    if the_vmag_spread gt vmag_spread then begin
        lprmsg, 'vmag spread is too large, skip ...'
        return, retval
    endif

    return, the_dir

end

;---Add test times.
    test_event_times = time_double([$
        '2007-11-20/17:18:10', $
        '2008-01-09/11:27:45', $
        '2008-01-19/12:02:55', $
        '2008-02-29/08:26:50', $
        '2014-08-28/10:10:40', $
        '2014-12-26/01:05:25', $
        '2016-10-13/12:22:35', $
        '2016-12-11/09:46:35', $
        '2017-03-28/03:00:40'])

    event_list = list()
    foreach test_event_time, test_event_times do event_list.add, test_event_time


;---Load all candidates.
    search_step = 'uniq_subgroup'
    candidates = azim_df_search_all_events(search_step=search_step, project=project)


;---Select the test events.
    index = list()
    foreach test_event_time, test_event_times do begin
        foreach candidate, candidates, id do begin
            if test_event_time eq candidate.time_range[0] then begin
                index.add, id
                lprmsg, 'Found the test event ...'
                break
            endif
        endforeach
    endforeach
    index = index.toarray()
    test_candidates = candidates[index]


;---Loop through.
    events = list()
    rejected = list()
    ;foreach candidate, test_candidates do begin
    foreach candidate, candidates do begin
        dir = azim_df_analyze_direction(candidate)
        if dir eq '' then begin
            candidate.direction = 'rejected'
            rejected.add, candidate
        endif else begin
            candidate.direction = dir
            events.add, candidate
        endelse
    endforeach
    stop


;---Copy files.
    root_dir = join_path([project.plot_dir,'diagnostic_plot'])
    in_dir = join_path([root_dir,'azim_df_uniq_df'])
    if file_test(in_dir) eq 0 then azim_df_subgroup_gen_diagnostic_plot, candidates, dirname='azim_df_uniq_df', project=project

    foreach event, events do begin
        out_dir = join_path([root_dir,'azim_df_analyze_direction',event.direction])
        if file_test(out_dir) eq 0 then file_mkdir, out_dir

        file_suffix = 'azim_df_event_'+strjoin(time_string(event.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
        in_file = join_path([in_dir,file_suffix])
        out_file = join_path([out_dir,file_suffix])
        file_copy, in_file, out_file, /overwrite
    endforeach

    out_dir = join_path([root_dir,'azim_df_analyze_direction','rejected'])
    if file_test(out_dir) eq 0 then file_mkdir, out_dir
    foreach event, rejected do begin
        file_suffix = 'azim_df_event_'+strjoin(time_string(event.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
        in_file = join_path([in_dir,file_suffix])
        out_file = join_path([out_dir,file_suffix])
        file_copy, in_file, out_file, /overwrite
    endforeach

end
