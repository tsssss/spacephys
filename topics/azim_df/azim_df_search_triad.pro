;+
; This is a newer version, which is more general than azim_df_search_candidate.
;
; Search for candidates with enough triads.
;
; Write settings to project, and write results to a text file.
;-

function azim_df_search_triad, search_setting, project=project, $
    reset=reset, test_time=test_time

;test_time = time_double('2014-08-28/10:00')
;test_time = time_double('2007-12-08/11:00')
    search_midn_candidate = project.search_midn_candidate
    retval = !null
    tab = constant('4space')

;---Settings for triad.
    the_key = 'search_triad'
    if ~search_setting.haskey(the_key) then message, 'No settings for triad ...'
    search_info = search_setting[the_key]

    if ~search_midn_candidate.haskey(the_key) then search_midn_candidate[the_key] = dictionary()
    settings = search_midn_candidate[the_key]
    small_angle_limit = 15.
    triad_angle_range = [small_angle_limit,180-small_angle_limit]
    default_settings = dictionary($
        'small_angle_limit', small_angle_limit, $
        'triad_angle_range', triad_angle_range, $
        'triad_min_duration', constant('secofhour'), $
        'file_suffix', project.name+'_search_midn_candidate_search_triad.txt')


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting triad search ...'
        file_delete, out_file, /allow_nonexistent
        lprmsg, 'Clear memory ...'
        del_data, '*'
    endif
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, roi_file
    endif

    if file_test(out_file) eq 0 then begin
        candidate_id = 0


        triad_angle_range = search_info.triad_angle_range
        triad_min_duration = search_info.triad_min_duration
        nvertex = 3
        lprmsg, 'Search candidates with good triads ...'
        lprmsg, tab+'Angle range (deg): ['+strjoin(string(triad_angle_range,format='(I0)'),',')+']'
        lprmsg, tab+'Min duration (min): '+string(triad_min_duration/60,format='(I0)')


    ;---Read ROI candidates, loop through each search type.
        roi_candidates = azim_df_search_roi(search_setting, project=project, reset=reset)
        probe_infos = project.probe_infos
        search_name = search_info.name
        lprmsg, ''
        lprmsg, 'Searching name: '+search_name+' ...'

        ; All probes to be searched.
        all_probes = search_setting.probes
        nall_probe = n_elements(all_probes)
        probe_labels = strarr(nall_probe)
        probe_colors = fltarr(nall_probe)
        foreach probe, all_probes, ii do begin
            probe_labels[ii] = strupcase(probe_infos[probe].short_name)
            probe_colors[ii] = probe_infos[probe].color
        endforeach
        lprmsg, tab+'Probes: '+strjoin(all_probes, ',')+' ...'

        ; Time range to be searched.
        full_time_range = search_setting.time_range
        time_step = project.search_candidate_time_step
        common_times = make_bins(full_time_range, time_step)
        ncommon_time = n_elements(common_times)
        lprmsg, tab+'Time range: '+strjoin(time_string(full_time_range,tformat='YYYY-MM-DD'), ' to ')+' ...'

        ; Min count for good triad.
        triad_min_count = search_info.triad_min_count
        lprmsg, tab+'Minimum # of good triad: '+string(triad_min_count,format='(I0)')


    ;---Load xxx_r_sm.
        lprmsg, tab+'Loading orbit data ...'
        ndim = 3
        r_sms = fltarr(ncommon_time,ndim,nall_probe)
        foreach probe, all_probes, ii do begin
            lprmsg, tab+tab+'Processing '+strupcase(probe)+' ...'
            azim_df_read_data, 'r_sm', time_range=full_time_range, probe=probe, project=project
            prefix = probe_infos[probe].prefix
            ; G14 is OK now, as the ROI time ranges already exclude when G14 do not have data.
            r_sms[*,*,ii] = get_var_data(prefix+'r_sm', at=common_times)
        endforeach
        r_sms[*,2,*] = 0    ; project to x-y plane.



    ;---Search for each region.
        region_name = search_setting.name
        lprmsg, ''
        lprmsg, tab+'Search in '+region_name+' ...'

    ;---Prepare data for calc triad.
        roi_flags = intarr(ncommon_time)
        probe_flags = intarr(ncommon_time,nall_probe)
        foreach candidate, roi_candidates do begin
            if candidate.search_name ne search_name then continue
            if candidate.region ne region_name then continue
            for section_id=0, candidate.nsection-1 do begin
                section_time_range = candidate.time_range_list[section_id]
                section_probes = candidate.probe_list[section_id]
                index = where_pro(common_times,'[]', section_time_range, count=count)
                if count eq 0 then message, 'Inconsistency, stop here ...'
                roi_flags[index] = 1
                foreach probe, section_probes do probe_flags[index,(where(all_probes eq probe))[0]] = 1
            endfor
        endforeach
        roi_index = where(roi_flags eq 1, nroi_time)
        if nroi_time eq 0 then message, 'Inconsistency, stop here ...'
        roi_times = common_times[roi_index]
        roi_flags = roi_flags[roi_index]
        probe_flags = probe_flags[roi_index,*]
        roi_rsms = r_sms[roi_index,*,*]

        roi_combos = choose_from(all_probes, nvertex)
        foreach combo, roi_combos do combo = combo[sort(combo)]
        if search_name eq 'beyond_15Re' then begin
            foreach combo, roi_combos, combo_id do if strjoin(combo,'_') eq strjoin('th'+['a','d','e'],'_') then roi_combos.remove, combo_id
        endif
        nroi_combo = n_elements(roi_combos)
        roi_triad_flags = intarr(nroi_time,nroi_combo)

    ;---Loop through time and calc triad.
        foreach candidate, roi_candidates do begin
            if candidate.search_name ne search_name then continue
            if candidate.region ne region_name then continue
            if keyword_set(test_time) then if product(candidate.time_range-test_time) lt 0 then stop

            for section_id=0, candidate.nsection-1 do begin
                section_time_range = candidate.time_range_list[section_id]
                section_probes = candidate.probe_list[section_id]
                time_index = where_pro(roi_times,'[]', section_time_range, count=ntime)
                if ntime eq 0 then message, 'Inconsistency, stop here ...'
                the_combos = (search_name eq 'beyond_15Re')? roi_combos: choose_from(section_probes, nvertex)
                foreach combo, the_combos do begin
                    combo_index = (roi_combos.where(combo))[0]
                    probe_index = intarr(nvertex)
                    foreach probe, combo, ii do probe_index[ii] = where(all_probes eq probe)
                    triad_angles = triangle_angles(roi_rsms[time_index,*,probe_index])
                    triad_flags = intarr(ntime)+1
                    for ii=0, nvertex-1 do begin
                        index = where_pro(triad_angles[*,ii], '][', triad_angle_range, count=count)
                        if count ne 0 then triad_flags[index] = 0
                    endfor
                    roi_triad_flags[time_index,combo_index] = triad_flags
                endforeach
            endfor
        endforeach


    ;---Summarize to get the candidates.
        roi_triad_counts = intarr(ncommon_time)
        roi_triad_counts[roi_index] = total(roi_triad_flags,2)
        roi_triad_count_var = search_name+'_'+region_name+'_triad_count'
        store_data, roi_triad_count_var, common_times, roi_triad_counts
        add_setting, roi_triad_count_var, /smart, {$
            display_type: 'scalar', $
            short_name: search_name+'!C  '+region_name+'!C  triad #', $
            unit: '(#)'}

        index = where(roi_triad_counts ge triad_min_count, count)
        if count eq 0 then message, 'No candidate, stop here ...'
        times = common_times[index]
        time_ranges = time_to_range(times, time_step=time_step)
        durations = time_ranges[*,1]-time_ranges[*,0]
        index = where(durations ge triad_min_duration, count)
        if count eq 0 then message, 'No candidate, stop here ...'
        time_ranges = time_ranges[index,*]


    ;---Break down into pieces when probes are fixed.
        lprmsg, tab+tab+'Writing results to file: '+out_file+' ...'
        ftouch, out_file
        ntime_range = n_elements(time_ranges)/2
        for ii=0, ntime_range-1 do begin
            time_range = reform(time_ranges[ii,*])
            if keyword_set(test_time) then if product(time_range-test_time) lt 0 then stop
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

            index = where_pro(roi_times,'[]',time_range)
            the_triad_flags = roi_triad_flags[index,*]
            the_times = roi_times[index]
            combo_counts = total(the_triad_flags, 1)
            index = where(combo_counts ne 0, navailable_combo)
            if navailable_combo lt triad_min_count then message, 'Inconsistency, stop here ...'
            the_triad_flags = the_triad_flags[*,index]
            the_combos = roi_combos[index]

            ; Divide into pieces according to combo.
            boundarys = time_range
            diff = fix(abs(the_triad_flags[1:-1,*]-the_triad_flags[0:-2,*]))
            diff = total(diff,2)
            index = where(diff ge 1, count)
            if count ne 0 then boundarys = [boundarys, the_times[index+1]]
            boundarys = sort_uniq(boundarys)
            nsection = n_elements(boundarys)-1
            for jj=0, nsection-1 do begin
                the_time_range = boundarys[jj:jj+1]
                the_flags = the_triad_flags[where_pro(the_times,'[)',the_time_range),*]
                available_combos = the_combos[where(total(the_flags,1) ne 0)]
                available_probes = list()
                foreach combo, available_combos do available_probes.add, combo, /extract
                available_probes = sort_uniq(available_probes.toarray())
                candidate.time_range_list.add, the_time_range
                candidate.probe_list.add, available_probes
            endfor

            ; Get all probes for all sections.
            all_section_probes = []
            for jj=0, nsection-1 do all_section_probes = [all_section_probes,candidate.probe_list[jj]]
            all_section_probes = sort_uniq(all_section_probes)
            candidate['all_probes'] = all_section_probes

            ; Combine as long as probe list are the same.
            probe_lists = strarr(nsection)
            for jj=0, nsection-1 do probe_lists[jj] = strjoin(sort_uniq(candidate.probe_list[jj]),'_')
            section_ids = intarr(nsection)
            section_id = 0
            for jj=1, nsection-1 do begin
                if probe_lists[jj] ne probe_lists[jj-1] then section_id += 1
                section_ids[jj] = section_id
            endfor
            nuniq_section = max(section_ids)+1
            time_range_list = list()
            probe_list = list()
            for jj=0, nuniq_section-1 do begin
                index = where(section_ids eq jj)
                the_time_range = []
                foreach the_index, index do the_time_range = [the_time_range,candidate.time_range_list[the_index]]
                probe_list.add, strsplit(probe_lists[index[0]],'_',/extract)
                time_range_list.add, minmax(the_time_range)
            endfor
            candidate.nsection = nuniq_section
            candidate.probe_list = probe_list
            candidate.time_range_list = time_range_list

            ; Write to file.
            azim_df_search_candidate_write_file, candidate, filename=out_file
        endfor

        lprmsg, tab+'Updating project ...'
        project.done_search_midn_candidate_search_triad = 1
        update_project, project

        lprmsg, tab+'Done search triad ...'
    endif


    lprmsg, ''
    lprmsg, 'Read triad candidate from file ...'
    candidates = azim_df_search_candidate_read_file(out_file)
    return, candidates

end