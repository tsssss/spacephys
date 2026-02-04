
pro micro_injection_spatial_distribution_manual_time_ranges_v01


    test = 0

    google_dir = googledir()
    data_dir = join_path([google_dir,'works','2024_microinjection','data'])
    plot_dir = join_path([google_dir,'works','2024_microinjection','plot'])
    base = 'micro_injection_survey_round1_event_list.txt'
    manual_tr_file = join_path([data_dir,base])
    lines = read_all_lines(manual_tr_file)
    nline = n_elements(lines)
    secofday = constant('secofday')

    trs = dblarr(nline,2)
    for ii=0,nline-1 do begin
        info = lines[ii]
        date = strmid(info,0,9)
        times = strmid(info,10,11)
        tr_str = date+'/'+strsplit(times,'-',extract=1)
        the_tr = time_double(tr_str,tformat='YYYY_MMDD/hh:mm')
        if the_tr[1] le the_tr[0] then begin
            the_tr[1] += secofday
        endif
        trs[ii,*] = the_tr
    endfor

    mission_probe = 'mms4'
    xrange = [1,-1]*15
    yrange = [1,-1]*15
    plot_file = join_path([plot_dir,'micro_injection_spatial_distribution_manual_time_ranges_v01.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, size=[8,8]
    tpos = sgcalcpos(1)
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, $
        ystyle=1, yrange=yrange, $
        position=tpos, nodata=1, iso=1
    
    tmp = lets_add_earth()
    
    psym = 6
    symsize = 0.1
    for ii=0,nline-1 do begin
        tr = reform(trs[ii,*])
        r_var = lets_read('orbit', tr, source=['mms','4'], coord='sm')
        r_sm = get_var_data(r_var, in=tr)
        xxs = r_sm[*,0]
        yys = r_sm[*,1]
        plots, xxs, yys, psym=psym, symsize=symsize
    endfor


    if keyword_set(test) then stop
    sgclose


end