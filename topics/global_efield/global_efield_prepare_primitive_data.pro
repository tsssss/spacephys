;+
; 1. Prepare year-based tplot files, for existing <mission>_load_<routine_name>
; 2. Combine the data to a cdf file.
;
; routine_name. A string to map to <mission>_read_<routine_name>, e.g., 'density', 'bfield', etc.
; probe=. A string to set mission_probe, e.g., 'polar', 'rbspa', 'thd', etc.
;-

;+
; Internal program.
;-
pro global_efield_prepare_primitive_data_internal, $
    routine_suffix, tplot_var_suffix, settings, $
    probe=mission_probe, project=project, errmsg=errmsg

;---Internal routine. Do not check inputs.
    project_updated = 0

;---Check if var exists in the combined_file
    probe_info = project[mission_probe]
    the_var = probe_info.prefix+tplot_var_suffix
    combined_file = join_path([project.data_dir,probe_info.combined_file_suffix])
    if cdf_has_var(the_var, filename=combined_file) then return

;---Check if the year-based files exist.
    the_key = routine_suffix+'_year_based_suffix'
    if ~probe_info.haskey(the_key) then begin
        years = probe_info['years']
        year_based_suffix = []
        foreach year, years do begin
            the_base = mission_probe+'_'+routine_suffix+'_'+string(year,format='(I4)')+'.tplot'
            the_suffix = join_path([routine_suffix,the_base])
            year_based_suffix = [year_based_suffix,the_suffix]
        endforeach
        probe_info[the_key] = year_based_suffix
        project_updated = 1
    endif
    if project_updated then update_project, project
    year_based_suffix = probe_info[the_key]
    year_based_files = year_based_suffix
    foreach file, year_based_files, ii do year_based_files[ii] = join_path([project.data_dir,file])


    ; Loop through each year-based file.
    years = probe_info['years']
    routine_name = probe_info.routine_name+'_read_'+routine_suffix
    probe = probe_info['probe']
    time_step = project['time_step']
    secofday = 86400d
    foreach the_file, year_based_files, ii do begin
        if file_test(the_file) then continue
        try_to_load_year_based_file = 1
        time_range = time_double(string(fix(years[ii])+[0,1],format='(I4)'))

    ;---Clear memory.
        del_data, the_var

    ;---Load data to memory through <mission>_read_<var>.
        if routine_suffix eq 'orbit' then begin
            call_procedure, routine_name, time_range, probe=probe
        endif else begin
            ndate = total(time_range*[-1,1])/secofday
            dates = smkarthm(time_range[0],secofday,ndate, 'x0')
            ptr_time = ptrarr(ndate)
            ptr_data = ptrarr(ndate)
            ndate_time = round(secofday/time_step)
        ;---Load data by day, for each day, coerce to the common time_step.
        ; To control the memory size.
            foreach date, dates, jj do begin
                date_time_range = date+[0,secofday]
                call_procedure, routine_name, date_time_range, probe=probe, errmsg=errmsg
                if errmsg ne '' then continue
                if tnames(the_var) eq '' then continue
                get_data, the_var, uts, dat
                if n_elements(uts) eq 1 then continue
                if n_elements(data_dims) eq 0 then data_dims = size(dat[0,*],/dimensions)
                date_time = smkarthm(date_time_range[0],time_step,ndate_time, 'x0')
                data_dims[0] = n_elements(date_time)
                date_data = make_array(dimension=data_dims, /float, value=!values.f_nan)
                for kk=0, ndate_time-1 do begin
                    index = lazy_where(uts, date_time[kk]+[0,time_step], count=count)
                    if count eq 0 then continue
                    tdata = dat[index,*]
                    index = where(finite(snorm(tdata)), count)  ; snorm is used to eliminate all dimension if one dimension is nan.
                    if count eq 0 then continue
                    vals = date_data[0,*]
                    foreach val, vals, ll do vals[ll] = mean(tdata[index,ll])
                    date_data[kk,*] = vals
                endfor
                ptr_time[jj] = ptr_new(date_time)
                ptr_data[jj] = ptr_new(date_data)
            endforeach

        ;---Combine data of all days into one file of the year.
            nrecs = lon64arr(ndate)
            for jj=0, ndate-1 do begin
                if ~ptr_valid(ptr_time[jj]) then continue
                nrecs[jj] = n_elements(*ptr_time[jj])
            endfor
            nrec = total(nrecs)
            if nrec eq 0 then continue
            time = dblarr(nrec)
            data_dims[0] = nrec
            data = make_array(dimension=data_dims, /float, value=!values.f_nan)

            i0 = 0
            for jj=0, ndate-1 do begin
                ; No data for the current day.
                if nrecs[jj] eq 0 then continue
                ; Fill the time and data in.
                i1 = i0+nrecs[jj]-1
                time[i0:i1] = *ptr_time[jj]
                data[i0:i1,*] = *ptr_data[jj]
                ; Free memory.
                *ptr_time[jj] = !null
                *ptr_data[jj] = !null
                i0 += nrecs[jj]
            endfor

            for jj=0, ndate-1 do begin
                if ~ptr_valid(ptr_time[jj]) then continue else ptr_free, ptr_time[jj]
                if ~ptr_valid(ptr_data[jj]) then continue else ptr_free, ptr_data[jj]
            endfor

            store_data, the_var, time, data
        endelse

    ;---Save data of the year from memory to disk.
        if tnames(the_var) eq '' then continue
        tplot_save, the_var, filename=the_file
    endforeach

    ; Remove the year_based_files do not exist.
    ; This happens when there is no data for a whole year.
    ; Remove the file suppress loading data for that year in the future.
    if keyword_set(try_to_load_year_based_file) then begin
        errmsg = ''
        flags = list()
        foreach the_file, year_based_files do flags.add, file_test(the_file)
        index = where(flags.toarray() eq 1, count)
        if count eq 0 then begin
            errmsg = handle_error('No data for the given var:'+the_var+' ...')
            probe_info.remove, the_key
        endif else begin
            probe_info[the_key] = year_based_suffix[index]
        endelse
        update_project, project
    endif


;---Combine data in year-based files and save to CDF.
    if ~probe_info.haskey(the_key) then return
    year_based_suffix = probe_info[the_key]
    year_based_files = year_based_suffix
    foreach file, year_based_files, ii do year_based_files[ii] = join_path([project.data_dir,file])
    time_var = 'ut_sec'

    time = []
    data = []
    foreach file, year_based_files do begin
        tplot_restore, filename=file
        get_data, the_var, uts, dat
        time = [time,uts]
        data = [data,dat]
    endforeach
    data = float(data)
    index = uniq(time)
    time = time[index]
    data = data[index,*,*,*,*,*,*,*]

    if routine_suffix eq 'orbit' then begin
        time_var_settings = dictionary($
            'unit', 'sec', $
            'time_var_type', 'unix')
        cdf_save_var, time_var, value=time, filename=combined_file
        cdf_save_setting, time_var_settings, filename=combined_file, varname=time_var
    endif else begin
        common_times = cdf_read_var(time_var, filename=combined_file)
        ncommon_time = n_elements(common_times)
        data_dims = size(data[0,*], /dimensions)
        data_dims[0] = ncommon_time
        common_data = make_array(dimension=data_dims, value=!values.f_nan)
        common_time_range = common_times[[0,ncommon_time-1]]
        index = lazy_where(time, '[]', common_time_range, count=count)
        time = time[index]
        data = data[index,*,*,*,*,*,*,*]
        time_index = round((time-common_time_range[0])/time_step)
        common_data[time_index,*,*,*,*,*,*,*] = data
        data = temporary(common_data)
    endelse

    dep_key = 'depend_0'
    if ~settings.haskey(dep_key) then settings[dep_key] = time_var
    cdf_save_var, the_var, value=data, filename=combined_file
    cdf_save_setting, settings, filename=combined_file, varname=the_var

end






;+
; The overarching program.
;-
pro global_efield_prepare_primitive_data, routine_name, probe=mission_probe, project=project

    if n_elements(routine_name) eq 0 then message, 'No routine_name ...'
    if n_elements(mission_probe) eq 0 then message, 'No mission_probe ...'
    if n_elements(project) eq 0 then project = global_efield_load_project()
    if ~project.haskey(mission_probe) then message, 'Invalid mission_probe ...'
    probe_info = project[mission_probe]
    project_updated = 0

;---Check if combined file exists.
    the_key = 'combined_file_suffix'
    if ~probe_info.haskey(the_key) then begin
        combined_file_suffix = join_path(['combined_data',project.name+'_'+mission_probe+'_combined_data_v01.cdf'])
        probe_info[the_key] = combined_file_suffix     ; dict is ref-passed.
        project_updated = 1
    endif
    if project_updated then update_project, project
    combined_file = join_path([project.data_dir,probe_info[the_key]])

    ; Create the combined file and write orbit data to it.
    if file_test(combined_file) eq 0 then begin
        settings = dictionary($
            'project', project.name, $
            'author', 'Sheng Tian', $
            'file_type', 'combined_data_file')
        cdf_save_setting, settings, filename=combined_file
        global_efield_prepare_primitive_data, 'orbit', probe=mission_probe, project=project
    endif


;---Dispatch.
    case routine_name of
        'orbit': begin
            tplot_var_suffix = 'r_gsm'
            routine_suffix = 'orbit'
            settings = dictionary($
                'display_type', 'vector', $
                'short_name', 'R', $
                'unit', 'Re', $
                'coord', 'GSM', $
                'coord_type', 'car', $
                'coord_labels', ['x','y','z'])
            end
        'density': begin
            tplot_var_suffix = 'ele_n'
            routine_suffix = 'density'
            settings = dictionary($
                'display_type', 'scalar', $
                'short_name', 'N!S!Uele!R!N', $
                'unit', 'cm!U-3!N', $
                'ylog', 1)
            end
        'ion_vel': begin
            tplot_var_suffix = 'u_gsm'
            routine_suffix = 'ion_vel'
            settings = dictionary($
                'display_type', 'vector', $
                'short_name', 'U!S!Uion!R!N', $
                'unit', 'km/s', $
                'coord', 'GSM', $
                'coord_type', 'car', $
                'coord_labels', ['x','y','z'])
            end
        'bfield': begin
            tplot_var_suffix = 'b_gsm'
            routine_suffix = 'bfield'
            settings = dictionary($
                'display_type', 'vector', $
                'short_name', 'B', $
                'unit', 'nT', $
                'coord', 'GSM', $
                'coord_type', 'car', $
                'coord_labels', ['x','y','z'])
            end
        'efield': begin
            tplot_var_suffix = 'e_gsm'
            routine_suffix = 'efield'
            settings = dictionary($
                'display_type', 'vector', $
                'short_name', 'E', $
                'unit', 'mV/m', $
                'coord', 'GSM', $
                'coord_type', 'car', $
                'coord_labels', ['x','y','z'])
        end
        'pflux': begin
            tplot_var_suffix = 'pf_gsm_norm'
            routine_suffix = 'pflux'
            settings = dictionary($
                'display_type', 'vector', $
                'short_name', 'S', $
                'unit', 'mW/m!U2!N', $
                'coord', 'GSM', $
                'coord_type', 'car', $
                'coord_labels', ['x','y','z'])
        end
        'ion_temp': begin
            tplot_var_suffix = 'ion_t'
            routine_suffix = 'ion_temp'
            settings = dictionary($
                'display_type', 'scalar', $
                'short_name', 'T!S!Uion!R!N', $
                'unit', 'eV', $
                'ylog', 1)
            end
        'ele_temp': begin
            tplot_var_suffix = 'ele_t'
            routine_suffix = 'ele_temp'
            settings = dictionary($
                'display_type', 'scalar', $
                'short_name', 'T!S!Uele!R!N', $
                'unit', 'eV', $
                'ylog', 1)
            end
    endcase

;---Do some real work.
    global_efield_prepare_primitive_data_internal, $
        routine_suffix, tplot_var_suffix, settings, $
        probe=mission_probe, project=project

;---Check and fix bad data.
    case routine_suffix of
        'orbit': global_efield_detect_and_fix_bad_orbit_data, mission_probe, project=project
        else: ; do nothing.
    endcase

end

project = global_efield_load_project()
mission_probes = project.all_mission_probes
mission_probes = 'c'+['1','2','3','4']
mission_probes = 'th'+letters('e')
vars = ['efield']

mission_probes = 'rbsp'+['b']
vars = 'pflux'
foreach mission_probe, mission_probes do foreach var, vars do global_efield_prepare_primitive_data, var, probe=mission_probe, project=project
end
