;+
; This is a newer version, which is more general than azim_df_search_candidate.
;
; Search for candidates in required ROI
; Write results to a text file and return.
;-

function azim_df_search_roi, search_setting, project=project, $
    reset=reset, test_time=test_time

;test_time = time_double('2014-08-28/10:00')
    retval = !null
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Settings for search DF in ROI.
    the_key = 'search_roi'
    if ~search_setting.haskey(the_key) then message, 'No setting for '+the_key+'...'
    search_info = search_setting[the_key]


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting '+the_key+'...'
        file_delete, out_file, /allow_nonexistent
    endif
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif

    if file_test(out_file) eq 0 then begin
        candidate_id = 0

        roi_min_duration = project.roi_min_duration
        overall_roi = project.overall_roi
        pdyn = overall_roi.pdyn
        mlt_range = search_info.mlt_range
        region_name = search_info.name

        lprmsg, 'Searching candidates in ROI ...'
        lprmsg, tab+'Region: '+region_name
        lprmsg, tab+'MLT range (hr): ['+strjoin(string(mlt_range,format='(I0)'),',')+']'
        lprmsg, tab+'Magnetopause model Pdyn (nPa) = '+string(pdyn,format='(I0)')
        lprmsg, tab+'Min duration (min): '+string(roi_min_duration/60,format='(I0)')
            ; Required # of probes in ROI.


    ;---Loop through each search type.
        probe_infos = project.probe_infos
        search_types = project.search_types
        foreach search_type, search_types do begin
            search_name = search_type.name

            ; All probes to be searched.
            probes = search_type.probes
            nprobe = n_elements(probes)
            probe_labels = strarr(nprobe)
            probe_colors = fltarr(nprobe)
            foreach probe, probes, ii do begin
                probe_labels[ii] = strupcase(probe_infos[probe].short_name)
                probe_colors[ii] = probe_infos[probe].color
            endforeach
            lprmsg, tab+'Probes: '+strjoin(probes, ',')+' ...'

            probe_min_count = search_type.roi_min_count

            ; Time range to be searched.
            full_time_range = search_type.time_range
            time_step = 300.
            common_times = make_bins(full_time_range, time_step)
            ncommon_time = n_elements(common_times)
            azim_df_load_basic_data, project=project
            lprmsg, tab+'Time range: '+strjoin(time_string(full_time_range,tformat='YYYY-MM-DD'), ' to ')+' ...'

        ;---Flags for each probe when it is in the "region of interest".
            foreach probe, probes do begin
                lprmsg, tab+tab+'Processing '+strupcase(probe)+' ...'
                prefix = probe_infos[probe].prefix
                roi_flags = bytarr(ncommon_time)+1

                r_sm = get_var_data(prefix+'r_sm', at=common_times)
                r_gsm = cotran(r_sm, common_times, 'sm2gsm')

                ; magnetopause.
                index = where(check_if_in_magn(r_gsm, dynamic_pressure=pdyn) eq 0, count)
                if count ne 0 then roi_flags[index] = 0

                ; rxy.
                rxy_range = overall_roi.rxy_range
                if search_name eq 'beyond_15Re' then begin
                    if probe eq 'thb' or probe eq 'thc' then rxy_range >= 15.
                endif
                rxy = snorm(r_sm[*,0:1])
                index = where_pro(rxy, '][', rxy_range, count=count)
                if count ne 0 then roi_flags[index] = 0

                ; "mlt" in SM, but really should be in MAG.
                mlt = get_var_data(prefix+'pseudo_mlt', at=common_times)
                index = where_pro(mlt, '][', mlt_range, count=count)
                if count ne 0 then roi_flags[index] = 0

                ; Stay in ROI longer than roi_min_duration.
                index = where(roi_flags eq 1, count)
                if count ne 0 then begin
                    time_ranges = time_to_range(common_times[index], time_step=time_step)
                    durations = time_ranges[*,1]-time_ranges[*,0]
                    index = where(durations ge roi_min_duration, count)
                    if count ne 0 then begin
                        roi_flags[*] = 0
                        time_ranges = time_ranges[index,*]
                        ntime_range = n_elements(time_ranges)/2
                        for ii=0, ntime_range-1 do roi_flags[where_pro(common_times,'[]',reform(time_ranges[ii,*]))] = 1
                    endif
                endif

                the_var = prefix+'roi_flag'
                store_data, the_var, common_times, roi_flags
                add_setting, the_var, /smart, {$
                    display_type: 'scalar', $
                    yrange: [-0.1,1.1], $
                    short_name: strupcase(probe_infos[probe].short_name), $
                    unit: '#'}
            endforeach

        ;---Combine all flags into one.
            roi_flag_var = 'roi_flag_'+search_name+'_'+region_name
            roi_flags = intarr(ncommon_time, nprobe)
            foreach probe, probes, ii do begin
                prefix = probe_infos[probe].prefix
                roi_flags[*,ii] = get_var_data(prefix+'roi_flag')
            endforeach
            store_data, roi_flag_var, common_times, roi_flags, limits={$
                ytitle:+region_name+'!CROI flag (#)', colors:probe_colors, labels:probe_labels, yrange:[-0.1,1.1]}

        ;---Count of # of probes in ROI at each time tag.
            roi_counts = total(roi_flags,2)
            roi_count_var = 'roi_count_'+search_name+'_'+region_name
            store_data, roi_count_var, common_times, roi_counts, limits={$
                ytitle:region_name+'!CROI count (#)', labels:region_name}


        ;---Find the times when there are enough probes, and all required probes are present.
            lprmsg, tab+tab+'Checking # of probes in ROI ...'
            time_index = roi_counts ge probe_min_count
            index = where(time_index eq 1, count)
            if count eq 0 then message, 'No candidate, stop here ...'
            if n_elements(required_probes) ne 0 then begin
                lprmsg, tab+tab+'Checking the existence of required probes ...'
                nrequired_probe = n_elements(required_probes)
                probe_index = intarr(nrequired_probe)
                foreach probe, required_probes, ii do probe_index[ii] = where(probes eq probe)
                roi_required_counts = total(roi_flags[*,probe_index],2)
                time_index = time_index and (roi_required_counts eq nrequired_probe)
                index = where(time_index eq 1, count)
                if count eq 0 then lprmsg, 'No candidate, stop here ...'
            endif
            times = common_times[index]


        ;---Rule out time ranges that are too short.
            lprmsg, tab+tab+'Excluding time ranges that are too short ...'
            time_ranges = time_to_range(times, time_step=time_step)
            durations = time_ranges[*,1]-time_ranges[*,0]
            index = where(durations gt roi_min_duration, count)
            if count eq 0 then lprmsg, 'No candidate, stop here ...'
            time_ranges = time_ranges[index,*]


        ;---Break down into sections with fixed probes.
            lprmsg, tab+tab+'Writing results to file: '+out_file+' ...'
            ftouch, out_file
            ntime_range = n_elements(time_ranges)/2
            for ii=0, ntime_range-1 do begin
                time_range = reform(time_ranges[ii,*])
                if keyword_set(test_time) then if product(candidate.time_range-test_time) lt 0 then stop
                candidate_id += 1
                candidate = dictionary($
                    'id', candidate_id, $
                    'time_range', time_range, $
                    'duration', total(time_range*[-1,1])/60, $  ; in min.
                    'nsection', 0, $
                    'search_name', search_name, $
                    'region', region_name, $
                    'all_probes', !null, $
                    'time_range_list', list(), $
                    'probe_list', list())

                index = where_pro(common_times,'[]',time_range)
                the_roi_flags = roi_flags[index,*]
                the_times = common_times[index]
                probe_counts = total(the_roi_flags, 1)
                index = where(probe_counts ne 0, navailable_probe)
                if navailable_probe lt probe_min_count then message, 'Inconsistency, stop here ...'
                the_roi_flags = the_roi_flags[*,index]

                boundarys = time_range
                if navailable_probe gt probe_min_count then begin
                    diff = fix(abs(the_roi_flags[1:-1,*]-the_roi_flags[0:-2,*]))
                    diff = total(diff,2)
                    index = where(diff ge 1, count)
                    if count ne 0 then boundarys = [boundarys, the_times[index+1]]
                endif
                boundarys = sort_uniq(boundarys)
                nsection = n_elements(boundarys)-1
                time_range_list = list()
                probe_list = list()
                for jj=0, nsection-1 do begin
                    the_time_range = boundarys[jj:jj+1]
                    the_roi_flags = roi_flags[where_pro(common_times,'[)',the_time_range),*]
                    probe_counts = total(the_roi_flags, 1)
                    available_probes = probes[where(total(the_roi_flags,1) ne 0)]
                    time_range_list.add, the_time_range
                    probe_list.add, available_probes
                endfor
                candidate.nsection = nsection
                candidate.probe_list = probe_list
                candidate.time_range_list = time_range_list

                ; Get all probes for all sections.
                all_section_probes = []
                for jj=0, nsection-1 do all_section_probes = [all_section_probes,candidate.probe_list[jj]]
                all_section_probes = sort_uniq(all_section_probes)
                candidate['all_probes'] = all_section_probes

                ; Write to file.
                azim_df_search_candidate_write_file, candidate, filename=out_file
            endfor
        endforeach
    endif


    lprmsg, ''
    lprmsg, 'Read ROI candidate from file ...'
    roi_candidates = azim_df_search_candidate_read_file(out_file)
    return, roi_candidates

end

    search_info = dictionary($
        'name', 'pre_midn', $
        'mlt_range', [-9,0], $
        'file_suffix', 'azim_df_search_pre_midn_search_roi.txt')
    search_setting = dictionary('search_roi', search_info)
    candidates = azim_df_search_roi(search_setting, project=project)
end
