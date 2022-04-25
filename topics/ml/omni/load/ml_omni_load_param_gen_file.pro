;+
; Read omni parameters: dst, ae, solar wind parameters.
;-

pro ml_omni_load_param_gen_file, input_date, filename=out_file, resolution=time_step, downsample_method=downsample_method, _extra=extra

test = 0

    if n_elements(input_date) eq 0 then begin
        errmsg = 'No input time ...'
        return
    endif

    if n_elements(time_step) eq 0 then time_step = ml_time_step()
    ; can be mean, median, log_mean.
    if n_elements(downsample_method) eq 0 then downsample_method = 'mean'

    year_str = time_string(input_date[0],tformat='YYYY')
    year = float(year_str)
    time_range = time_double([year_str,string(year+1,format='(I04)')])
    secofday = constant('secofday')
if keyword_set(test) then time_range = time_range[0]+[0,30*secofday]
    days = make_bins(time_range+[0,-1]*secofday,secofday)
    common_times = make_bins(time_range+[0,-1]*time_step,time_step)+time_step*0.5

    prefix = 'omni_'

;---Init file.
    if file_test(out_file) eq 0 then cdf_touch, out_file

;---Add time.
    time_var = 'unix_time'
    if ~cdf_has_var(time_var, filename=out_file) then begin
        common_times = make_bins(time_range+[0,-1]*time_step, time_step)+time_step*0.5

        cdf_save_var, time_var, value=common_times, filename=out_file
        settings = dictionary($
            'FIELDNAM', 'Time', $
            'UNITS', 'sec', $
            'TEXT', 'center time of the sampling interval', $
            'VAR_TYPE', 'support_data' )
        cdf_save_setting, settings, var=time_var, filename=out_file
    endif
    common_times = cdf_read_var(time_var, filename=out_file)
    ntime = n_elements(common_times)
    fillval = !values.f_nan
    


;---Add omni parameters.
    omni_vars = prefix+['symh','symd','asyd','asyh','al','au',$
        'sw_'+['time_shift','b_gse','v_gse','p','n','t']]
    load_data = 0
    foreach var, omni_vars do begin
        if cdf_has_var(var, filename=out_file) then continue
        load_data = 1
        break
    endforeach
    if load_data then begin
        files = omni_load(time_range, id='cdaweb%hro2%1min')

        var_list = list()
        ; Vxyz is in GSE.
        ; AE = AU-AL, AU is often small, so AE and AL are similar.
        xyz = constant('xyz')
        in_vars = ['SYM_H','SYM_D','ASY_D','ASY_H',$
            'AL_INDEX','AU_INDEX','Timeshift',$
            ['BX','BY','BZ']+'_GSE',['Vx','Vy','Vz'],$
            'Pressure','proton_density','T']
        out_vars = prefix+['symh','symd','asyd','asyh',$
            'al','au','sw_time_shift',$
            'sw_b'+xyz+'_gse','sw_v'+xyz+'_gse',$
            'sw_p','sw_n','sw_t']
        var_list.add, dictionary($
            'in_vars', in_vars, $
            'out_vars', out_vars, $
            'time_var_name', 'Epoch', $
            'time_var_type', 'epoch' )
        read_vars, time_range, files=files, var_list=var_list
        foreach var, out_vars, var_id do begin
            vatt = cdf_read_setting(in_vars[var_id], filename=files[0])
            if vatt.haskey('FILLVAL') then begin
                val = vatt['FILLVAL']
                get_data, var, times, data
                data = float(data)
                index = where(data eq val, count)
                if count ne 0 then data[index] = fillval
                ; Need to rewrite b/c original data are in long.
                store_data, var, times, data
            endif
        endforeach
        
        var = prefix+'time_shift'
        get_data, var, times, data
        ;index = where(abs(data) gt 1e4, count)  ; valid max=1e4 from skeleton.
        index = where(finite(get_var_data(prefix+'sw_vx_gse'),nan=1), count)
        if count ne 0 then begin
            data[index] = fillval
            store_data, var, times, data
        endif
        

    ;---Snap data to common_times.
        ndim = 3

        ; Solar wind B GSE.
        the_var = prefix+'sw_b_gse'
        settings = dictionary($
            'FIELDNAM', 'Solar wind B GSE', $
            'UNITS', 'nT', $
            'VAR_TYPE', 'data', $
            'DEPEND_0', time_var )
        data = fltarr([ntime,ndim])+fillval
        get_data, prefix+'sw_bx_gse', times, xx
        get_data, prefix+'sw_by_gse', times, yy
        get_data, prefix+'sw_bz_gse', times, zz
        foreach time, common_times, time_id do begin
            time_index = where(times ge time-0.5*time_step and $
                    times le time+0.5*time_step, count)
            if count eq 0 then continue
            if downsample_method eq 'mean' then begin
                data[time_id,0] = mean(xx[time_index], nan=1)
                data[time_id,1] = mean(yy[time_index], nan=1)
                data[time_id,2] = mean(zz[time_index], nan=1)
            endif else if downsample_method eq 'median' then begin
                data[time_id,0] = median(xx[time_index], nan=1)
                data[time_id,1] = median(yy[time_index], nan=1)
                data[time_id,2] = median(zz[time_index], nan=1)
            endif
        endforeach
        cdf_save_var, the_var, value=data, filename=out_file
        cdf_save_setting, settings, var=the_var, filename=out_file

        ; Solar wind V GSE.
        the_var = prefix+'sw_v_gse'
        settings = dictionary($
            'FIELDNAM', 'Solar wind V GSE', $
            'UNITS', 'km/s', $
            'VAR_TYPE', 'data', $
            'DEPEND_0', time_var )
        data = fltarr([ntime,ndim])+fillval
        get_data, prefix+'sw_vx_gse', times, xx
        get_data, prefix+'sw_vy_gse', times, yy
        get_data, prefix+'sw_vz_gse', times, zz
        foreach time, common_times, time_id do begin
            time_index = where(times ge time-0.5*time_step and $
                    times le time+0.5*time_step, count)
            if count eq 0 then continue
            if downsample_method eq 'mean' then begin
                data[time_id,0] = mean(xx[time_index], nan=1)
                data[time_id,1] = mean(yy[time_index], nan=1)
                data[time_id,2] = mean(zz[time_index], nan=1)
            endif else if downsample_method eq 'median' then begin
                data[time_id,0] = median(xx[time_index], nan=1)
                data[time_id,1] = median(yy[time_index], nan=1)
                data[time_id,2] = median(zz[time_index], nan=1)
            endif
        endforeach
        cdf_save_var, the_var, value=data, filename=out_file
        cdf_save_setting, settings, var=the_var, filename=out_file

        ; Solar wind pressure.
        ; Solar wind proton density.
        ; Solar wind temperature (what species??).
        ; Solar wind time shift.
        vars = prefix+'sw_'+['p','n','t','time_shift']
        labels = 'Solar wind '+['dynamic pressure','proton density','temperature','time shift']
        units = ['nPa','cm!U-3','eV','sec']
        foreach the_var, vars, var_id do begin
            settings = dictionary($
                'FIELDNAM', labels[var_id], $
                'UNITS', units[var_id], $
                'VAR_TYPE', 'data', $
                'DEPEND_0', time_var )
            data = fltarr(ntime)+fillval
            get_data, the_var, times, xx
            if the_var eq prefix+'sw_t' then xx *= (1d/11600)   ; convert K to eV.
            foreach time, common_times, time_id do begin
                time_index = where(times ge time-0.5*time_step and $
                        times le time+0.5*time_step, count)
                if count eq 0 then continue
                if downsample_method eq 'mean' then begin
                    data[time_id] = mean(xx[time_index], nan=1)
                endif else if downsample_method eq 'median' then begin
                    data[time_id] = median(xx[time_index], nan=1)
                endif
            endforeach
            cdf_save_var, the_var, value=data, filename=out_file
            cdf_save_setting, settings, var=the_var, filename=out_file
        endforeach

        ; Sym-H, Sym-D, Asy-H, Asy-D, AU, AL.
        vars = prefix+['symh','symd','asyh','asyd','au','al']
        labels = ['Sym-H','Sym-D','Asy-H','Asy-D','AU','AL']+' index'
        foreach the_var, vars, var_id do begin 
            settings = dictionary($
                'FIELDNAM', labels[var_id], $
                'UNITS', 'nT', $
                'VAR_TYPE', 'data', $
                'DEPEND_0', time_var )
            data = fltarr(ntime)+fillval
            get_data, the_var, times, xx
            foreach time, common_times, time_id do begin
                time_index = where(times ge time-0.5*time_step and $
                        times le time+0.5*time_step, count)
                if count eq 0 then continue
                if downsample_method eq 'mean' then begin
                    data[time_id] = mean(xx[time_index], nan=1)
                endif else if downsample_method eq 'median' then begin
                    data[time_id] = median(xx[time_index], nan=1)
                endif
            endforeach
            cdf_save_var, the_var, value=data, filename=out_file
            cdf_save_setting, settings, var=the_var, filename=out_file
        endforeach
    endif



;---Add omni hourly data.
    omni_vars = prefix+['f107','pc_index']
    load_data = 0
    foreach var, omni_vars do begin
        if cdf_has_var(var, filename=out_file) then continue
        load_data = 1
        break
    endforeach
    if load_data then begin
        files = omni_load(time_range, id='cdaweb%hourly')

        var_list = list()
        in_vars = ['PC_N_INDEX','F10_INDEX']
        out_vars = prefix+['pc_index','f107']
        var_list.add, dictionary($
            'in_vars', in_vars, $
            'out_vars', out_vars, $
            'time_var_name', 'Epoch', $
            'time_var_type', 'epoch' )

        read_vars, time_range, files=files, var_list=var_list
        foreach var, out_vars, var_id do begin
            vatt = cdf_read_setting(in_vars[var_id], filename=files[0])
            if vatt.haskey('FILLVAL') then begin
                val = vatt['FILLVAL']
                get_data, var, times, data
                data = float(data)
                index = where(data eq val, count)
                if count ne 0 then data[index] = fillval
                ; Need to rewrite b/c original data are in long.
                store_data, var, times, data
            endif
        endforeach
        

        the_var = prefix+'f107'
        settings = dictionary($
            'FIELDNAM', 'F10.7', $
            'UNITS', '10!U-22!N J/s-m!U2!N-Hz', $
            'VAR_TYPE', 'data', $
            'DEPEND_0', time_var )
        get_data, the_var, times, data
        index = where(times mod secofday eq 0, count)
        if count ne 0 then begin
            times = times[index]+secofday*0.5
            data = data[index]
            data = interpol(data, times, common_times, nan=1)
        endif else begin
            data = fltarr(ntime)+fillval
        endelse
        cdf_save_var, the_var, value=data, filename=out_file
        cdf_save_setting, settings, var=the_var, filename=out_file

        the_var = prefix+'pc_index'
        settings = dictionary($
            'FIELDNAM', 'F10.7', $
            'UNITS', '#', $
            'VAR_TYPE', 'data', $
            'DEPEND_0', time_var )
        get_data, the_var, times, data
        data = interpol(data, times, common_times, nan=1)
        cdf_save_var, the_var, value=data, filename=out_file
        cdf_save_setting, settings, var=the_var, filename=out_file
    endif


end

stop
time_range = ['1996','2020']
resolutions = [60d,300d]
foreach resolution, resolutions do begin
    files = ml_omni_load_param(time_range, resolution=resolution)
    foreach file, files do begin
        base = file_basename(file)
        date = strmid(base,16,4)
        ml_omni_load_param_gen_file, date, filename=file
    endforeach
endforeach

end