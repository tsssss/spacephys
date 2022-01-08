;+
; Check if a given candidate is a coherent DF.
;-


function azim_df_filter_coherent_df, candidate, project=project, $
    settings=settings, $
    log_file=log_file

    retval = dictionary()
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(settings) eq 0 then begin
        settings = dictionary($
            'max_angle_diff', 15., $
            'max_vmag_ratio', 0.2)
    endif


    max_angle_diff = settings.max_angle_diff
    max_vmag_ratio = settings.max_vmag_ratio
    if n_elements(candidate) eq 0 then begin
        lprmsg, 'Settings for filtering triad ...', log_file
        lprmsg, 'Max angle diff (deg):'+string(max_angle_diff,format='(I0)'), log_file
        lprmsg, 'Max vmag diff ratio (#): '+string(max_vmag_ratio,format='(F3.1)'), log_file
        lprmsg, '', log_file
        return, settings
    endif

    retval = dictionary()
    if n_elements(project) eq 0 then project = azim_df_load_project()
;    azim_df_load_basic_data, project=project, scale_width=scale_width
    lprmsg, '', log_file
    msg = 'Processing candidate: '+strjoin(time_string(candidate.time_range),' to ')
    lprmsg, msg, log_file


;---Check vmag.
    triad_list = candidate.triad_list
    ntriad = triad_list.length
    vmags = fltarr(ntriad)
    foreach triad, triad_list, ii do vmags[ii] = triad.vmag_obs_time
    lprmsg, tab+'vmags (km/s): '+strjoin(string(vmags,format='(I0)'),','), log_file
    vmag_mean = mean(vmags)
    vmag_stddev = stddev(vmags)
    vmag_ratio = vmag_stddev/vmag_mean
    vmag_ratio = round(vmag_ratio*10)/10.
    lprmsg, tab+'vmag_ratio = vmag_stddev/mag_mean: '+$
        string(vmag_ratio,format='(F4.1)')+' = '+$
        string(vmag_stddev,format='(I0)')+'/'+$
        string(vmag_mean,format='(I0)'), log_file
    if vmag_ratio gt max_vmag_ratio then begin
        lprmsg, 'Inconsistent vmag, skip ...', log_file
        return, retval
    endif

;---Check vhat.
    vhat_angles = fltarr(ntriad)
    foreach triad, triad_list, ii do vhat_angles[ii] = atan(triad.vhat_obs_time[1],triad.vhat_obs_time[0])
    vhat_angles *= constant('deg')
    lprmsg, tab+'vhat angles (deg): '+strjoin(string(vhat_angles,format='(I0)'),','), log_file
    vhat_angle_stddev = stddev(vhat_angles)
    vhat_angle_stddev = round(vhat_angle_stddev)
    lprmsg, tab+'vhat_angle_stddev (deg): '+string(vhat_angle_stddev,format='(I0)'), log_file
    if vhat_angle_stddev gt max_angle_diff then begin
        lprmsg, 'Inconsistent vhat, skip ...', log_file
        return, retval
    endif

    msg = 'Pass ...'
    lprmsg, msg, log_file
    return, candidate

end
