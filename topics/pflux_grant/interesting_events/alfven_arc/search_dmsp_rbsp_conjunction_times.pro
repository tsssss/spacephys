;+
; Search DMSP and RBSP conjunctions.
;-

function search_dmsp_rbsp_conjunction_times_read_info, file

    lines = read_all_lines(file)
    nline = n_elements(lines)
    nsection = 16
    nevent = float(nline)/nsection
    infos = orderedhash()
    for ii=0,nevent-1 do begin
        i0 = ii*nsection
        i1 = i0+nsection-1
        the_lines = lines[i0:i1]
        j0 = 1
        id = the_lines[j0+0]
        time_range = time_double(strsplit(the_lines[j0+1],' ',extract=1),tformat='YYYY_MMDD_hhmm')
        rbsp_probe = (strsplit(the_lines[j0+2],':',extract=1))[1]
        rbsp_hem = (strsplit(the_lines[j0+3],':',extract=1))[1]
        rbsp_mlt_if = (strsplit(the_lines[j0+4],':',extract=1))[1]
        rbsp_mlat_if = (strsplit(the_lines[j0+5],':',extract=1))[1]
        rbsp_mlt_range = (strsplit(the_lines[j0+6],':',extract=1))[1]
        rbsp_mlat_range = (strsplit(the_lines[j0+7],':',extract=1))[1]
        rbsp_dis = (strsplit(the_lines[j0+8],':',extract=1))[1]

        dmsp_probe = (strsplit(the_lines[j0+9],':',extract=1))[1]
        dmsp_hem = (strsplit(the_lines[j0+10],':',extract=1))[1]
        dmsp_mlt_if = (strsplit(the_lines[j0+11],':',extract=1))[1]
        dmsp_mlat_if = (strsplit(the_lines[j0+12],':',extract=1))[1]
        dmsp_mlt_range = (strsplit(the_lines[j0+13],':',extract=1))[1]
        dmsp_mlat_range = (strsplit(the_lines[j0+14],':',extract=1))[1]
        

        info = dictionary($
            'time_range', time_range, $
            'id', id )
        info['rbsp'] = dictionary($
            'probe', rbsp_probe, $
            'mlat_range', float(strsplit(rbsp_mlat_range,',',extract=1)), $
            'mlat_if', float(strsplit(rbsp_mlat_if,',',extract=1)), $
            'mlt_range', float(strsplit(rbsp_mlt_range,',',extract=1)), $
            'mlt_if', float(strsplit(rbsp_mlt_if,',',extract=1)), $
            'dis', float(rbsp_dis), $
            'hemisphere', rbsp_hem )
        info['dmsp'] = dictionary($
            'probe', dmsp_probe, $
            'mlat_range', float(strsplit(dmsp_mlat_range,',',extract=1)), $
            'mlat_if', float(strsplit(dmsp_mlat_if,',',extract=1)), $
            'mlt_range', float(strsplit(dmsp_mlt_range,',',extract=1)), $
            'mlt_if', float(strsplit(dmsp_mlt_if,',',extract=1)), $
            'hemisphere', dmsp_hem )
        infos[id] = info
    endfor
    return, infos

end

function search_dmsp_rbsp_conjunction_times_write_info, info, filename=file

    msgs = list()
    msgs.add, ''
    msgs.add, info.id
    msgs.add, strjoin(time_string(info.time_range,tformat='YYYY_MMDD_hhmm'),' ')

    msgs.add, 'RBSP:'+info.rbsp.probe
    msgs.add, 'hemisphere:'+info.rbsp.hemisphere
    msgs.add, 'mlt_i,f(h):'+strjoin(string(info.rbsp.mlt_if,format='(F5.1)'),',')
    msgs.add, 'mlat_i,f(deg):'+strjoin(string(info.rbsp.mlat_if,format='(F5.1)'),',')
    msgs.add, 'mlt_range(h):'+strjoin(string(info.rbsp.mlt_range,format='(F5.1)'),',')
    msgs.add, 'mlat_range(deg):'+strjoin(string(info.rbsp.mlat_range,format='(F5.1)'),',')
    msgs.add, 'dis(Re):'+string(info.rbsp.dis,format='(F3.1)')

    msgs.add, 'DMSP:'+info.dmsp.probe
    msgs.add, 'hemisphere:'+info.dmsp.hemisphere
    msgs.add, 'mlt_i,f(h):'+strjoin(string(info.dmsp.mlt_if,format='(F5.1)'),',')
    msgs.add, 'mlat_i,f(deg):'+strjoin(string(info.dmsp.mlat_if,format='(F5.1)'),',')
    msgs.add, 'mlt_range(h):'+strjoin(string(info.dmsp.mlt_range,format='(F5.1)'),',')
    msgs.add, 'mlat_range(deg):'+strjoin(string(info.dmsp.mlat_range,format='(F5.1)'),',')


    foreach msg, msgs do lprmsg, msg, file
    return, msgs

end

function search_dmsp_rbsp_conjunction_times, input_time_range, $
    filename=file, dmlt=dmlt, $
    rbsp_probe=rbsp_probe, dmsp_probe=dmsp_probe, $
    errmsg=errmsg

    errmsg = ''
    retval = !null

    if n_elements(input_time_range) eq 0 then input_time_range = ['2013','2015']
    if n_elements(rbsp_probe) eq 0 then rbsp_probe = 'a'
    if n_elements(dmsp_probe) eq 0 then dmsp_probe = 'f18'
    time_range = time_double(input_time_range)
    if n_elements(dmlt) eq 0 then dmlt = 1.


;---RBSP location.
    ; Load data.
    rbsp_pad_time = 4.5*3600
    rbsp_time_range = time_range+[-1,1]*rbsp_pad_time
    rbsp_r_var = rbsp_read_orbit(rbsp_time_range, probe=rbsp_probe, get_name=1)
    if check_if_update(rbsp_r_var, rbsp_time_range) then rbsp_r_var = rbsp_read_orbit(rbsp_time_range, probe=rbsp_probe)
    rbsp_mlat_var = rbsp_read_mlat(rbsp_time_range, probe=rbsp_probe, get_name=1)
    if check_if_update(rbsp_mlat_var, rbsp_time_range) then rbsp_mlat_var = rbsp_read_mlat(rbsp_time_range, probe=rbsp_probe)
    rbsp_mlt_var = rbsp_read_mlt(rbsp_time_range, probe=rbsp_probe, get_name=1)
    if check_if_update(rbsp_mlt_var, rbsp_time_range) then rbsp_mlt_var = rbsp_read_mlt(rbsp_time_range, probe=rbsp_probe)
    get_data, rbsp_mlt_var, times, mlts
    index = where(mlts ge 12, count)
    if count ne 0 then begin
        mlts[index] -= 24
        store_data, rbsp_mlt_var, times, mlts
    endif
    
    ; Get the time of interest: around apogee.
    diss = snorm(get_var_data(rbsp_r_var, times=times))
    min_dis = 4
    index = where(diss ge min_dis, count)
    if count eq 0 then return, retval
    rbsp_time_ranges = times[time_to_range(index,time_step=1)]

    ; Remove partial orbit.
    perigee_dis = 2.
    index = where(diss le perigee_dis, count)
    if count eq 0 then message, 'Inconsistency ...'
    perigee_time_ranges = times[time_to_range(index,time_step=1)]
    perigee_times = total(perigee_time_ranges,2)*0.5
    index = where(rbsp_time_ranges[*,0] gt min(perigee_times) $
        and rbsp_time_ranges[*,1] lt max(perigee_times))
    rbsp_time_ranges = rbsp_time_ranges[index,*]
    nrbsp_time_range = n_elements(rbsp_time_ranges[*,0])
    

;---DMSP location.
    infos = orderedhash()
    abs_mlat_range = [60,70]
    for rbsp_id=0,nrbsp_time_range-1 do begin
        dmsp_time_range = reform(rbsp_time_ranges[rbsp_id,*])
        dmsp_mlat_vars = dmsp_read_mlat_vars(dmsp_time_range, probe=dmsp_probe, errmsg=errmsg)
        if errmsg ne '' then continue
        dmsp_mlat_var = dmsp_mlat_vars[0]
        mlats = get_var_data(dmsp_mlat_var, times=times)
        index = where_pro(abs(mlats), '[]', abs_mlat_range, count=count)
        if count eq 0 then continue
        dmsp_time_ranges = times[time_to_range(index,time_step=1)]
        ndmsp_time_range = n_elements(dmsp_time_ranges[*,0])
        dmsp_mlt_var = dmsp_mlat_vars[1]
        get_data, dmsp_mlt_var, times, mlts
        index = where(mlts ge 12, count)
        if count ne 0 then begin
            mlts[index] -= 24
            store_data, dmsp_mlt_var, times, mlts
        endif
        for dmsp_id=0,ndmsp_time_range-1 do begin
            the_tr = reform(dmsp_time_ranges[dmsp_id,*])
            the_duration = abs(total(the_tr*[-1,1]))
            if the_duration le 60 then continue
            ; round to minutes.
            if (the_tr[1] mod 60) eq 0 then the_tr[1]-= 60
            the_tr = the_tr-(the_tr mod 60)+[0,60]
            ; get dmsp and rbsp locations.
            dmsp_mlts = get_var_data(dmsp_mlt_var, in=the_tr)
            dmsp_mlats = get_var_data(dmsp_mlat_var, in=the_tr)
            dmsp_mlt = mean(dmsp_mlts)
            ndmsp = n_elements(dmsp_mlts)
            rbsp_mlts = get_var_data(rbsp_mlt_var, in=the_tr)
            rbsp_mlt = mean(rbsp_mlts)
            nrbsp = n_elements(rbsp_mlts)

            print, 'DMSP '+strupcase(dmsp_probe)
            print, 'Time range (UT): '+strjoin(time_string(the_tr),',')
            print, 'MLat range (deg): '+strjoin(string(minmax(dmsp_mlats),format='(F5.1)'),',')
            print, 'MLT range (h): '+strjoin(string(minmax(dmsp_mlts),format='(F5.1)'),',')
            
            the_dmlt = abs(dmsp_mlt-rbsp_mlt)
            if the_dmlt gt dmlt then continue
            
            rbsp_mlats = get_var_data(rbsp_mlat_var, in=the_tr)
            dmsp_hem = (mean(dmsp_mlats) ge 0)? 'north': 'south'
            rbsp_hem = (mean(rbsp_mlats) ge 0)? 'north': 'south'
            info = dictionary('time_range', the_tr)
            info['rbsp'] = dictionary($
                'probe', rbsp_probe, $
                'mlat_range', minmax(rbsp_mlats), $
                'mlat_if', rbsp_mlats[[0,nrbsp-1]], $
                'mlt_range', minmax(rbsp_mlts), $
                'mlt_if', rbsp_mlts[[0,nrbsp-1]], $
                'dis', snorm(get_var_data(rbsp_r_var, at=mean(the_tr))), $
                'hemisphere', rbsp_hem )
            info['dmsp'] = dictionary($
                'probe', dmsp_probe, $
                'mlat_range', minmax(dmsp_mlats), $
                'mlat_if', dmsp_mlats[[0,ndmsp-1]], $
                'mlt_range', minmax(dmsp_mlts), $
                'mlt_if', dmsp_mlts[[0,ndmsp-1]], $
                'hemisphere', dmsp_hem )
            info['id'] = strjoin([info.rbsp.probe,info.dmsp.probe,$
                time_string(info.time_range,tformat='YYYY_MMDD_hhmm')], '%')

            msgs = search_dmsp_rbsp_conjunction_times_write_info(info, filename=file)
            infos[info.id] = info
        endfor
    endfor

    return, infos

end



;---Settings.
    data_dir = join_path([srootdir(),'data'])
    years = make_bins([2013,2014],1)
    rbsp_probes = ['a','b']
    dmsp_probes = 'f'+['16','17','18']

;---Step1: generate conjunction for each year.
    step1_files = list()
    foreach year, years do begin
        search_time_range = time_double(string(year+[0,1],format='(I4)'))
        foreach rbsp_probe, rbsp_probes do begin
            foreach dmsp_probe, dmsp_probes do begin
                base = 'search_'+string(year,format='(I4)')+$
                    '_rbsp'+rbsp_probe+'_dmsp'+dmsp_probe+'_v02.txt'
                file = join_path([data_dir,base])
                step1_files.add, file
                if file_test(file) eq 1 then continue
                ftouch, file
                print, file
                ;continue
                infos = search_dmsp_rbsp_conjunction_times(search_time_range, $
                    rbsp_probe=rbsp_probe, dmsp_probe=dmsp_probe, filename=file)
            endforeach
        endforeach
    endforeach

;---Step2: parse step1 files and filter down.
    mlt_range = [-4,4]
    ;mlt_range = [-1,1]*4
    min_dst = -40

    base = 'search_dmsp_rbsp_mlt_range_'+strjoin(string(mlt_range,format='(I0)'),'_')+'_dst_'+string(min_dst,format='(I0)')+'_v01.txt'
    step2_file = join_path([data_dir,base])
    if file_test(step2_file) then stop


    infos = orderedhash()
    foreach file, step1_files do begin
        the_infos = search_dmsp_rbsp_conjunction_times_read_info(file)
        foreach info, the_infos do begin
            if infos.haskey(info.id) then continue
            infos[info.id] = info
        endforeach
    endforeach
    ninfo = infos.count()
    print, ninfo
    
    

    mlt_flags = fltarr(ninfo)
    foreach id, infos.keys(), ii do begin
        info = infos[id]
        
        ; MLT range.
        rbsp_mlt = mean(info.rbsp.mlt_range)
        index = where_pro(rbsp_mlt, '[]', mlt_range, count=count)
        if count eq 0 then begin
            mlt_flags[ii] = 1
            continue
        endif
        
;        ; DMSP orbit, want cross instead of skim
;        dmsp_mlat_range = minmax(info.dmsp.mlat_range)
;        dmsp_mlat_if = minmax(info.dmsp.mlat_if)
;        dmlat = abs(total(dmsp_mlat_range-dmsp_mlat_if))
;        if dmlat ge 0.1 then begin
;            mlt_flags[ii] = 1
;            continue
;        endif
        
        ; Same hemisphere.
        dmsp_hem = info.dmsp.hemisphere
        rbsp_hem = info.rbsp.hemisphere
        if dmsp_hem ne rbsp_hem then begin
            mlt_flags[ii] = 1
            continue
        endif
    endforeach
    index = where(mlt_flags eq 1, count)
    if count ne 0 then begin
        infos.remove, (infos.keys())[index]
    endif
    ninfo = infos.count()
    
    
    full_time_range = time_double(string(minmax(years)+[0,1],format='(I4)'))+[-1,1]*constant('secofday')
    dst_var = 'dst'
    if check_if_update(dst_var, full_time_range) then begin
        omni_read_index, full_time_range
        options, dst_var, 'requested_time_range', full_time_range
    endif
    dst_flags = fltarr(ninfo)
    foreach id, infos.keys(), ii do begin
        info = infos[id]
        dst_flags[ii] = mean(get_var_data(dst_var, in=info.time_range+[-1,1]*600))
    endforeach
    index = where(dst_flags gt min_dst, count)
    if count ne 0 then begin
        infos.remove, (infos.keys())[index]
    endif
    ninfo = infos.count()
    print, ninfo
    

    ftouch, step2_file
    foreach info, infos do begin
        tmp = search_dmsp_rbsp_conjunction_times_write_info(info, filename=step2_file)
    endforeach

stop
;tr = ['2013-05-01','2013-05-02']
;infos = search_dmsp_rbsp_conjunction_times(tr, rbsp_probe='b', dmsp_probe='f18')


end