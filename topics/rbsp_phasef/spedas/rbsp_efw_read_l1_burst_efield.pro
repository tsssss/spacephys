;+
; Read EFW burst data. Save data in rbspx_efw_eb1_<coord>.
;
; datatype can be 'vb1','vb2'.
;
; coord=. A lower case string sets the coordinate of the output data. 'mgse' by default.
; keep_spin_axis=. A boolean, set to keep the spin axis E field. 0 by default.
;-

pro rbsp_efw_read_l1_burst_efield, tr, probe=probe, $
    datatype=datatype, trange=trange, $
    coord=coord, keep_spin_axis=keep_spin_axis, $
    level=level, verbose=verbose, downloadonly=downloadonly, $
    cdf_data=cdf_data,get_support_data=get_support_data, $
    tplotnames=tns, make_multi_tplotvar=make_multi_tplotvar, $
    varformat=varformat, valid_names=valid_names, files=files, $
    type=type, _extra=_extra

    rbsp_efw_init
    vb = keyword_set(verbose) ? verbose : 0
    vb = vb > !rbsp_efw.verbose

    if n_elements(coord) eq 0 then coord = 'mgse'
    if n_elements(keep_spin_axis) eq 0 then keep_spin_axis = 0
    if n_elements(datatype) eq 0 then datatype = 'vb1-split'

    data_types = ['vb1','vb2','vb1-split']
    data_type = datatype[0]
    index = where(data_types eq data_type, count)
    if count eq 0 then begin
        dprint, 'Invalid datatype: '+data_type+' ...', verbose=vb
        return
    endif
    if data_type eq 'vb1-split' then begin
        rbsp_efw_read_l1_burst_split, tr, probe=probe, $
            datatype=data_type, trange=trange, $
            level=level, verbose=verbose, downloadonly=downloadonly, $
            cdf_data=cdf_data,get_support_data=get_support_data, $
            tplotnames=tns, make_multi_tplotvar=make_multi_tplotvar, $
            varformat=varformat, valid_names=valid_names, files=files, $
            type=type, _extra=_extra
    endif else begin
        rbsp_efw_read_l1, tr, probe=probe, $
            datatype=data_type, trange=trange, $
            level=level, verbose=verbose, downloadonly=downloadonly, $
            cdf_data=cdf_data,get_support_data=get_support_data, $
            tplotnames=tns, make_multi_tplotvar=make_multi_tplotvar, $
            varformat=varformat, valid_names=valid_names, files=files, $
            type=type, _extra=_extra
    endelse


;---Calbirate E field.
    prefix = 'rbsp'+probe+'_'
    if data_type eq 'vb2' then begin
        data_type2 = 'vb2'
    endif else begin
        data_type2 = 'vb1'
    endelse
    v_var = prefix+'efw_'+data_type2
    get_data, v_var, data=dd
    if size(dd,/type) ne 8 then begin
        errmsg = 'No data ...'
        return
    endif
    store_data, v_var, dlimits={data_att:{units:'ADC'}}
    time_range = time_double(tr)
    timespan, time_range[0], total(time_range*[-1,1]), second=1
    rbsp_efw_cal_waveform, probe=probe, datatype=data_type2, trange=time_range

    ; Convert vsvy to esvy.
    cp0 = rbsp_efw_get_cal_params(time_range[0])
    cp = (probe eq 'a')? cp0.a: cp0.b
    boom_length = cp.boom_length
;    boom_shorting_factor = cp.boom_shorting_factor

    get_data, v_var, times, vsvy
    ntime = n_elements(times)
    ndim = 3
    esvy = dblarr(ntime,ndim)
    for eid=0,ndim-1 do begin
        vid = eid*2
        coef = 1d3/boom_length[eid]
        ;coef = 1d
        esvy[*,eid] = (vsvy[*,vid]-vsvy[*,vid+1])*coef
    endfor

    data_type3 = 'e'+strmid(data_type2,1,2)
    in_var = prefix+'efw_'+data_type3
    store_data, in_var, times, esvy

;---Remove DC offset, use E UVW.
    survey_time_range = minmax(times)+[-1,1]*30
    rbsp_efw_phasef_read_dc_offset, survey_time_range, probe=probe
    get_data, prefix+'efw_e_uvw_dc_offset', uts, euvw_offset
    euvw_offset[*,2] = 0    ; Do nothing to Ew.
    esvy -= sinterpol(euvw_offset, uts, times)
    store_data, in_var, times, esvy


;---Spin-axis data.
    if ~keep_spin_axis then esvy[*,2] = 0
    store_data, in_var, times, esvy


;---Convert vector from UVW to wanted coord.
    rgb = [6,4,2]
    xyz = ['x','y','z']
    get_data, in_var, times, vec, limits=lim
    vec = cotran(vec, times, 'uvw2'+coord[0], probe=probe, use_orig_quaternion=1)
    out_var = in_var+'_'+coord[0]
    store_data, out_var, times, vec, limits={$
        ytitle:prefix+'burst_efield!C[mV/m]', $
        labels:strupcase(coord)+' E'+xyz, $
        colors:rgb }

end


; Set the time and probe for loading data.
time_range = ['2013-06-10/05:57:20','2013-06-10/05:59:40']
probe = 'b'

time_range = ['2013-06-07/02:17:40','2013-06-07/02:19:20']
probe = 'a'

; Load the spinfit data.
rbsp_efw_read_l1_burst_efield, time_range, probe=probe, $
    datatype='vb1', coord='mgse'

prefix = 'rbsp'+probe+'_efw_'
vars = prefix+[$

    ; The E field in UVW.
    'eb1_mgse' ]

; Plot the variables.
tplot, vars, trange=time_range
end
