;+
; To study how to detect a random dipolarization.
;
; Use derivative of theta, then use gauss fit to find the arrival time, width and height.
;-

pro azim_df_study_detect_df, project=project, probe=probe, time_range=time_range

    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(probe) eq 0 then message, 'No input probe ...'
    if n_elements(time_range) ne 2 then message, 'Invalid input time_range ...'

;---Load the theta for given probe and time_range.
    prefix = probe+'_'
    search_name = (mean(time_range) le time_double('2010-01-01'))? 'beyond_15Re': 'within_15Re'
    search_settings = project.search_settings
    foreach search_setting, search_settings do if search_setting.name eq search_name then break
    data_file = join_path([project.data_dir,search_setting.data_file_suffix])

    theta_var = prefix+'theta'
    if check_if_update(theta_var, time_range) then azim_df_read_data_theta, probe=probe, time_range=time_range, filename=data_file
    b_gsm_var = prefix+'b_gsm'
    if check_if_update(b_gsm_var, time_range) then azim_df_read_data_b_gsm, probe=probe, time_range=time_range, filename=data_file


;---Preprocess theta.
    test_theta_var = prefix+'test_theta'
    get_data, theta_var, times, theta, limits=lim
    store_data, test_theta_var, times, theta, limits=lim
    ntime = n_elements(times)


    ; Rmove points when |B| is too small.
    ; Do not need this if we do the smooth.
    bmag_limit = 2. ; nT.
    bmag = snorm(get_var_data(b_gsm_var, times=times))
    index = where(bmag ge bmag_limit, count)
    if count lt ntime then begin
        theta = interpol(theta[index], times[index], times)
        store_data, test_theta_var, times, theta
    endif

    ; Smooth by boxcar mean at 80% ratio.
    time_step = project.time_step
    boxcar_window = 120.    ; sec.
    boxcar_width = boxcar_window/time_step
    boxcar_ratio = 0.8

    boxcar_boundary_times = make_bins(time_range, boxcar_window)
    nboxcar = n_elements(boxcar_boundary_times)-1
    boxcar_center_times = boxcar_boundary_times[0:nboxcar-1]+boxcar_window*0.5
    theta = get_var_data(test_theta_var)
    test_theta = fltarr(nboxcar)+!values.f_nan
    theta_stddev = fltarr(nboxcar)+!values.f_nan
    theta_upper = fltarr(nboxcar)+!values.f_nan
    theta_lower = fltarr(nboxcar)+!values.f_nan
    for ii=0, nboxcar-1 do begin
        index = where_pro(times, '[]', boxcar_boundary_times[ii:ii+1], count=count)
        if count eq 0 then continue
        the_theta = theta[index]
        index = where(finite(the_theta), count)
        if count lt 2 then continue
        the_theta = the_theta[index]
        index = sort(abs(the_theta))
        the_theta = the_theta[index]
        test_theta[ii] = mean(the_theta[0:count*boxcar_ratio])
        theta_stddev[ii] = stddev(the_theta[0:count*boxcar_ratio])
        theta_upper[ii] = max(the_theta[0:count*boxcar_ratio])
        theta_lower[ii] = min(the_theta[0:count*boxcar_ratio])
    endfor
;    test_theta = interpol(test_theta, boxcar_center_times, times, /quadratic)
    store_data, test_theta_var, boxcar_center_times, test_theta

    ; Standard deviation.
    min_theta = stddev(theta_stddev)
    theta_stddev = smooth(theta_stddev > min_theta, 7, /edge_truncate, /nan)

    smooth_theta = interpol(test_theta, boxcar_center_times, times, /quadratic)
    theta_stddev = interpol(theta_stddev, boxcar_center_times, times, /quadratic)
    theta_upper = interpol(theta_upper, boxcar_center_times, times, /quadratic)
    theta_lower = interpol(theta_lower, boxcar_center_times, times, /quadratic)
    theta_combo_var = prefix+'theta_combo'
    store_data, theta_combo_var, times, [[theta],[smooth_theta],[theta_upper],[theta_lower]], $
        limits={ytitle:'(deg)',labels:tex2str('theta')+['orig',' smooth',' upper', ' lower'],colors:sgcolor(['silver','red','tan','tan'])}


;---Calculate the derivative.
    dtheta_var = prefix+'dtheta'
    used_theta_var = test_theta_var
    ;used_theta_var = theta_var
    theta = get_var_data(used_theta_var, times=the_times)
    dtheta = deriv(the_times, theta)
    dtheta = interpol(dtheta, the_times, times, /quadratic)
    store_data, dtheta_var, times, dtheta
    add_setting, dtheta_var, /smart, {$
        display_type: 'scalar', $
        unit: 'deg/sec', $
        short_name: 'd'+tex2str('theta')+'/dt', $
        constant: 0}

;---Calcualte the local standard deviation of the derivative.
    dtheta_stddev = fltarr(nboxcar)+!values.f_nan
    for ii=0, nboxcar-1 do begin
        index = where_pro(times, '[]', boxcar_boundary_times[ii:ii+1], count=count)
        if count eq 0 then continue
        the_dtheta = dtheta[index]
        dtheta_stddev[ii] = stddev(the_dtheta, /nan)
    endfor
    nsigma = 3
    min_dtheta = stddev(dtheta_stddev)
    dtheta_stddev = smooth(dtheta_stddev > min_dtheta, 7, /edge_truncate, /nan)
    dtheta_stddev = interpol(dtheta_stddev, boxcar_center_times, times)

    dtheta_combo_var = prefix+'dtheta_combo'
    store_data, dtheta_combo_var, times, [[dtheta],[dtheta_stddev*nsigma],[dtheta_stddev]], $
        limits={ytitle:'(deg/sec)',labels:'d'+tex2str('theta')+['',' stddev x'+sgnum2str(nsigma),' stddev'],colors:sgcolor(['silver','red','blue']),constant:0}

test = 0
;---Select peaks in derivative that exceeds standard deviation.
    df_infos = list()
    index = where(dtheta gt dtheta_stddev, count)
    gauss_nterm = 3
    gauss_color = sgcolor('red')
    fwhm_coef = sqrt(-alog(0.5)*2)
    dtheta_step = stddev(deriv(dtheta))

    if count ne 0 then begin
        the_times = times[index]
        time_ranges = time_to_range(the_times, time_step=time_step)
        ntime_range = n_elements(time_ranges)/2

        if keyword_set(test) then begin
            window, /free, xsize=800, ysize=300
            plot, times-times[0], dtheta
        endif

        for ii=0, ntime_range-1 do begin
            ; Exclude the time ranges that are too short.
            the_time_range = reform(time_ranges[ii,*])
            index = where_pro(times, '[]', the_time_range, count=count)
            the_index = where(dtheta[index] ge dtheta_stddev[index]*nsigma, count)
            if count le gauss_nterm then continue


            ; Exclude the time ranges that are on the edges.
            if index[0] eq 0 then continue
            if index[-1] eq ntime-1 then continue

            yys = dtheta[index]
            xxs = times[index]
            yfit = gaussfit(xxs, yys, coeff, nterms=gauss_nterm, chisq=chisq)
            print, chisq
            if keyword_set(test) then begin
                oplot, xxs-times[0], yfit, color=gauss_color
                stop
            endif

;            arrival_time = coeff[1]
;            width = coeff[2]*fwhm_coef*2
;            ;if width/total(the_time_range*[-1,1]) ge 1.5 then continue
;            the_time_range = arrival_time+[-1,1]*width*0.5

            arrival_time = mean(the_time_range)
            the_theta = get_var_data(used_theta_var, at=the_time_range)
            the_value_range = minmax(the_theta)
            if product(the_value_range) gt 0 and min(abs(the_value_range)) ge min(interpol(theta_stddev,times,the_time_range)) then continue
            height = total(the_value_range*[-1,1])
            width = total(the_time_range*[-1,1])


            the_info = dictionary($
                'time_range', the_time_range, $
                'value_range', the_value_range, $
                'width', width, $
                'height', height, $
                'arrival_time', arrival_time)
            df_infos.add, the_info
        endfor
    endif

    yys = get_var_data(theta_var, in=time_range)
    yrange = minmax(yys)
    options, theta_combo_var, 'yrange', yrange
    options, theta_combo_var, 'ystyle', 1
    options, theta_combo_var, 'constant', 0

    sgopen, 0, xsize=8, ysize=8, /inch
    margins = [10,4,10,2]

    vars = prefix+['theta_combo','dtheta_combo','b_gsm']
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, margins=margins)
    tplot, vars, position=poss, trange=time_range

    xrange = time_range
    tpos = poss[*,0]
    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        /nodata, /noerase
    yys = get_var_data(theta_combo_var, in=time_range, times=xxs)
    yys = yys[*,1]
    bar_thick = 4
    foreach df_info, df_infos do begin
        the_time_range = df_info.time_range
        index = where_pro(xxs, '[]', the_time_range)
        the_xxs = xxs[index]
        the_yys = yys[index]
        plots, the_xxs, the_yys, color=sgcolor('blue'), thick=bar_thick
    endforeach
stop



end


if n_elements(project) eq 0 then project = azim_df_load_project()


    ; Fig 1-1 in Runov+2010.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-02-27/07:30','2009-02-27/09:00'])

;    ; Fig 1-2 in Runov+2010.
;    probes = 'th'+letters('e')
;    time_range = time_double(['2009-03-05/02:30','2009-03-05/04:00'])
;
;    ; Fig 1-3 in Runov+2010.
;    probes = 'th'+letters('e')
;    time_range = time_double(['2009-03-09/08:30','2009-03-09/10:00'])
;
;    ; Fig 2-1 in Runov+2010.
;    probes = 'th'+letters('e')
;    time_range = time_double(['2009-03-15/08:30','2009-03-15/09:30'])
;
;    ; Fig 2-2 in Runov+2010.
;    probes = 'th'+letters('e')
;    time_range = time_double(['2009-03-19/08:00','2009-03-19/09:00'])
;
;    ; Fig 2-3 in Runov+2010.
;    probes = 'th'+letters('e')
;    time_range = time_double(['2009-03-31/08:00','2009-03-31/09:00'])
;
;
;    ; The themis event.
;    time_range = time_double(['2008-02-29/08:00','2008-02-29/09:30'])
;    probes = 'th'+letters('e')

;    ; The Runov 2009 event.
;    time_range = time_double(['2009-02-27/07:30','2009-02-27/09:00'])
;    probes = 'th'+['b','c','d','e','a']

    ; My event.
    ;time_range = time_double(['2014-08-28/09:30','2014-08-28/11:00'])
    ;probes = ['rbspb','g13','g15','tha','thd','the']

;    time_range = time_double(['2016-10-13/12:00','2016-10-13/13:30'])
;    probes = ['rbspa','rbspb','g15','thd','g14','g13']




foreach probe, probes do azim_df_study_detect_df, project=project, probe=probe, time_range=time_range

vars = []
foreach probe, probes do vars = [vars,probe+'_'+['theta_combo','dtheta_combo']]
tplot, vars
end
