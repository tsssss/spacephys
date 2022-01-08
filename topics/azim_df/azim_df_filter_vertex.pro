;+
; Filter DFs with wanted shapes.
; Adopted from azim_df_filter_large_df.
;-

function azim_df_filter_vertex, df, project=project, settings=settings, $
    log_file=log_file

    tab = constant('4space')
    if n_elements(log_file) eq 0 then log_file = -1
    if n_elements(settings) eq 0 then begin
        settings = dictionary($
            'height_range', [0.,80], $
            'width_range', [1.,25]*60., $
            'min_scaled_height', 8., $
            'max_abs_scaled_theta', 80., $
            'scale_width', 4., $
            'section_time_range', [0.,10]*60., $
            'section_min_scaled_height', 0., $
            'section_min_ratio', 0.4, $
            'section_total_ratio', 0.8 )
    endif

    height_range = settings.height_range
    width_range = settings.width_range
    min_scaled_height = settings.min_scaled_height
    max_abs_scaled_theta = settings.max_abs_scaled_theta
    scale_width = settings.scale_width
    section_time_range = settings.section_time_range
    section_min_scaled_height = settings.section_min_scaled_height
    section_min_ratio = settings.section_min_ratio
    section_total_ratio = settings.section_total_ratio
    if n_elements(df) eq 0 then begin
        lprmsg, 'Settings for filtering dipolarizations ...', log_file
        lprmsg, tab+'Height range (deg): '+strjoin(string(height_range,format='(I0)'), ' to '), log_file
        lprmsg, tab+'Width range (min): '+strjoin(string(width_range/60,format='(I0)'), ' to '), log_file
        lprmsg, tab+'Minimum scaled height (deg): '+string(min_scaled_height,format='(I0)'), log_file
        lprmsg, tab+'Maximum abs scaled theta (deg): '+string(max_abs_scaled_theta,format='(I0)'), log_file
        lprmsg, tab+'Scale width: '+string(scale_width,format='(F3.1)'), log_file
        lprmsg, tab+'Section time range (sec): ['+strjoin(string(section_time_range,format='(I0)'),',')+']', log_file
        lprmsg, tab+'Section min scaled_height (deg): '+string(section_min_scaled_height,format='(I0)'), log_file
        lprmsg, tab+'Section total ratio (#): '+string(section_total_ratio,format='(F4.2)'), log_file
        lprmsg, tab+'Section continuous ratio (#): '+string(section_min_ratio,format='(F4.2)'), log_file
        lprmsg, '', log_file
        return, settings
    endif

    retval = dictionary()
    if n_elements(project) eq 0 then project = azim_df_load_project()
    azim_df_load_basic_data, project=project, scale_width=scale_width
    msg = 'Processing DF: '+string(df.probe,format='(A8)')+tab+time_string(df.obs_time)
    lprmsg, msg, log_file


;---Check width.
    msg = tab+'width (sec): '+string(df.width,format='(I0)')
    lprmsg, msg, log_file
    index = lazy_where(df.width, '][', width_range, count=count)
    if count ne 0 then begin
        msg = 'DF out of width range, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif

;---Check height.
    msg = tab+'height (deg): '+string(df.height,format='(F5.1)')
    lprmsg, msg, log_file
    index = lazy_where(df.height, '][', height_range, count=count)
    if count ne 0 then begin
        msg = 'DF out of height range, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif

;---Check scaled_height.
    df.scaled_height = azim_df_scale_theta(df.height, df.obs_mlt, width=scale_width)
    scaled_height = float(round(df.scaled_height))
    msg = tab+'scaled_height (deg): '+string(df.scaled_height,format='(F7.1)')
    lprmsg, msg, log_file
    if scaled_height le min_scaled_height then begin
        msg = 'DF out of scaled_height range, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif

    scaling_factor = df.scaled_height/df.height
    scaled_theta_ranges = df.theta_range*scaling_factor
    if max(abs(scaled_theta_ranges)) ge max_abs_scaled_theta then begin
        msg = 'Scaled theta range out of range, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif


;---Check shape.
    prefix = df.probe+'_'
    theta_var = prefix+'scaled_theta'
    the_time_range = df.obs_time+section_time_range
    scaled_heights = get_var_data(theta_var, in=the_time_range, times=times)
    index = where(scaled_heights ge section_min_scaled_height, count)
    if count eq 0 then begin
        lprmsg, 'Inconsistency in scaled_height ...', log_file
        return, retval
    endif
    common_data_rate = project.time_step
    time_ranges = time_to_range(times[index], time_step=common_data_rate)
    durations = time_ranges[*,1]-time_ranges[*,0]
    section_duration = total(section_time_range*[-1,1])

    ; Total duration.
    total_duration = total(durations)
    section_ratio = total_duration/section_duration
    section_ratio = round(section_ratio*10)/10.
    msg = tab+'section total ratio (#): '+string(section_ratio,format='(F3.1)')+$
        ', '+string(total_duration,format='(I5)')+' sec'
    lprmsg, msg, log_file
    if section_ratio ge section_total_ratio then begin
        msg = 'Pass ...'
        lprmsg, msg, log_file
        return, df
    endif

    ; Max duration.
    max_duration = max(durations, section_index)
    max_duration_time_range = reform(time_ranges[section_index,*])
    section_ratio = max_duration/section_duration
    section_ratio = round(section_ratio*10)/10.
    msg = tab+'section max ratio (#): '+string(section_ratio,format='(F3.1)')+$
        ', '+string(max_duration,format='(I5)')+$
        ' sec, '+strjoin(time_string(max_duration_time_range),' to ')
    lprmsg, msg, log_file
    if section_ratio lt section_min_ratio then begin
        msg = 'DF shape is bad, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif


    msg = 'Pass ...'
    lprmsg, msg, log_file
    return, df


;    if isa(df,'list') then begin
;        large_df = list()
;        foreach tmp, df do begin
;            tmp = azim_df_filter_vertex(tmp, project=project, $
;                settings=settings, log_file=log_file)
;            if n_elements(tmp) eq 0 then continue
;            large_df.add, tmp
;        endforeach
;        return, large_df
;    endif


end






    check_good = 1
    if check_good then begin
        df_strs = [$
            'thd       2008-02-29/08:27:00      480      23.97          24.2    -0.54    11.41    2008-02-29/08:23:00 2008-02-29/08:31:00     -13.96   10.01    -11.30,  1.61,  0.34', $

            'tha       2017-03-28/03:00:40      300      13.84          15.4    -1.87     8.75    2017-03-28/02:58:10 2017-03-28/03:03:10      -1.49   12.35     -7.72,  4.12,  0.59', $
            'the       2017-03-28/03:05:05      290      24.32          35.4    -3.47    11.35    2017-03-28/03:02:40 2017-03-28/03:07:30      -8.45   15.87     -6.98,  8.95,  1.03', $
            'thd       2017-03-28/03:08:45      130       6.26          10.2    -3.95    10.26    2017-03-28/03:07:40 2017-03-28/03:09:50      -1.20    5.06     -5.24,  8.82,  1.02', $
            'rbspa     2017-03-28/03:09:25      370       5.48          10.8    -4.66     5.73    2017-03-28/03:06:20 2017-03-28/03:12:30      -2.80    2.68     -1.97,  5.38, -0.19', $
            'g15       2017-03-28/03:14:05      670      11.49          33.8    -5.88     6.59    2017-03-28/03:08:30 2017-03-28/03:19:40      -8.76    2.73     -0.20,  6.59,  0.51', $

            'rbspa     2016-10-13/12:22:35      510      16.96          17.8     1.27     5.70    2016-10-13/12:18:20 2016-10-13/12:26:50      -8.96    8.00     -5.39, -1.86,  0.89', $
            'rbspb     2016-10-13/12:25:40      420      12.91          14.2     1.73     5.40    2016-10-13/12:22:10 2016-10-13/12:29:10      -6.46    6.45     -4.86, -2.36,  0.94', $
            'g15       2016-10-13/12:28:55      590      20.93          32.5     3.75     6.60    2016-10-13/12:24:00 2016-10-13/12:33:50      -9.74   11.19     -3.67, -5.48,  0.53', $
            'thd       2016-10-13/12:31:25      410      23.52          51.2     4.99     7.76    2016-10-13/12:28:00 2016-10-13/12:34:50     -15.29    8.23     -2.03, -7.49,  0.67', $
            'g14       2016-10-13/12:37:20      260       9.73          28.8     5.89     6.54    2016-10-13/12:35:10 2016-10-13/12:39:30      -7.74    1.98     -0.18, -6.54,  0.93', $
            'g13       2016-10-13/12:45:10      360      15.82         118.1     8.02     6.52    2016-10-13/12:42:10 2016-10-13/12:48:10      -6.83    8.99      3.29, -5.63,  1.11', $

            'thd       2009-04-07/07:05:15      670      33.50          34.6    -1.02    11.43    2009-04-07/06:59:40 2009-04-07/07:10:50     -13.36   20.14    -11.03,  3.01,  0.52', $
            'tha       2009-04-07/07:05:20      280      14.88          15.9    -1.43    11.67    2009-04-07/07:03:00 2009-04-07/07:07:40      -5.82    9.06    -10.87,  4.26, -0.72', $
            'thb       2009-04-07/07:06:40      660      52.31          67.5    -2.86    28.26    2009-04-07/07:01:10 2009-04-07/07:12:10     -26.15   26.16    -20.71, 19.23,  3.95', $
            'thc       2009-04-07/07:08:15      650      31.57          44.2    -3.28    18.57    2009-04-07/07:02:50 2009-04-07/07:13:40     -20.36   11.21    -12.13, 14.06,  1.32', $

            'thc       2008-01-09/11:27:45      330       8.14           8.7     1.49    19.27    2008-01-09/11:25:00 2008-01-09/11:30:30      -3.54    4.59    -17.83, -7.31,  0.73', $
            'tha       2008-01-09/11:32:35      390      16.45          17.9     1.62     9.85    2008-01-09/11:29:20 2008-01-09/11:35:50     -13.27    3.17     -8.98, -4.05,  0.22', $
            'the       2008-01-09/11:37:10      560      27.26          32.7     2.41    11.73    2008-01-09/11:32:30 2008-01-09/11:41:50     -16.21   11.05     -9.48, -6.91,  0.33', $
            'thd       2008-01-09/11:39:50      540      25.38          32.0     2.72    11.39    2008-01-09/11:35:20 2008-01-09/11:44:20     -12.78   12.59     -8.62, -7.45,  0.41', $

            'thd       2014-08-28/10:10:40      380      32.39          32.5     0.38     7.17    2014-08-28/10:07:30 2014-08-28/10:13:50     -21.68   10.71     -7.13, -0.71,  0.53', $
            'the       2014-08-28/10:13:20      880      28.80          30.1     1.19     9.32    2014-08-28/10:06:00 2014-08-28/10:20:40     -14.29   14.51     -8.87, -2.86,  1.00', $
            'g15       2014-08-28/10:15:15      450      31.60          32.8     1.11     6.59    2014-08-28/10:11:30 2014-08-28/10:19:00     -20.12   11.48     -6.31, -1.89,  0.56', $
            'tha       2014-08-28/10:20:25      350      10.71          11.9     1.85    11.12    2014-08-28/10:17:30 2014-08-28/10:23:20      -3.21    7.50     -9.84, -5.18,  1.56', $
            'rbspb     2014-08-28/10:44:55      370       3.41           8.6     5.45     5.46    2014-08-28/10:41:50 2014-08-28/10:48:00      -0.78    2.63     -0.78, -5.40,  1.90', $
            'g13       2014-08-28/10:47:25      430       6.35          17.0     5.62     6.51    2014-08-28/10:43:50 2014-08-28/10:51:00      -2.39    3.95     -0.64, -6.48,  1.11', $

            'thc       2008-03-20/12:06:00      340       9.31          12.8    -3.18    18.94    2008-03-20/12:03:10 2008-03-20/12:08:50      -2.95    6.37    -12.74, 14.01, -3.20', $
            'thb       2008-03-20/12:09:10      220       9.20          11.2    -2.50    28.13    2008-03-20/12:07:20 2008-03-20/12:11:00      -2.41    6.79    -22.32, 17.12, -3.36', $
            'thd       2008-03-20/12:16:35      470      26.34          26.4    -0.31     8.07    2008-03-20/12:12:40 2008-03-20/12:20:30     -16.90    9.44     -8.04,  0.66, -1.29', $
            'the       2008-03-20/12:17:30      500      19.10          19.5    -0.85     9.44    2008-03-20/12:13:20 2008-03-20/12:21:40      -7.29   11.81     -9.21,  2.08, -1.73']
        df_list = list()
        foreach df_str, df_strs do df_list.add, azim_df_vertex_read(df_str)
    endif else begin
        df_strs = [$
            'thc       2007-12-17/09:05:50      440      63.56         114.4     4.34    15.99    2007-12-17/09:02:10 2007-12-17/09:09:30       1.36   64.92     -6.75,-14.50,  2.35', $
            'tha       2014-08-28/10:42:50      180       7.50           8.4     1.91    11.30    2014-08-28/10:41:20 2014-08-28/10:44:20       1.21    8.71     -9.91, -5.42,  1.51', $
            'thb       2008-02-16/01:22:00      460      39.00          46.7    -2.40    19.50    2008-02-16/01:18:10 2008-02-16/01:25:50     -10.50   28.50    -15.79, 11.45,  3.37']
        df_list = list()
        foreach df_str, df_strs do df_list.add, azim_df_vertex_read(df_str)
    endelse


    log_file = -1
    foreach df, df_list do begin
        large_df = azim_df_filter_vertex(df, project=project, log_file=log_file)
        if check_good then begin
            if n_elements(large_df) eq 0 then stop
        endif else begin
            if n_elements(large_df) eq 1 then stop
        endelse
    endforeach
end
