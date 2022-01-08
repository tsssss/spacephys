;+
; Load data that have been coerced to the common time to memory.
;-

pro global_efield_prepare_bmod_t89_gsm, probe=mission_probe, project=project, filename=cdf_file

    model = 't89'
    coord_name = 'gsm'

    prefix = project[mission_probe].prefix
    var_suffix = 'r_gsm'
    global_efield_load_data, var_suffix, probe=mission_probe
    rgsm_var = prefix+var_suffix
    get_data, rgsm_var, times, rgsm, limits=lim

    coord_labels = ['x','y','z']
    ndim = n_elements(coord_labels)
    ntime = n_elements(times)
    b0gsm = fltarr(ntime,ndim)  ; in-situ B field in GSM.
    par = 2.
    for ii=0, ntime-1 do begin
        tilt = geopack_recalc(times[ii])
        ; in-situ position
        rx = rgsm[ii,0]
        ry = rgsm[ii,1]
        rz = rgsm[ii,2]
        ; in-situ B field.
        geopack_igrf_gsm, rx,ry,rz, bx,by,bz
        geopack_t89, par, rx,ry,rz, dbx,dby,dbz
        b0gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
    endfor
    data = b0gsm
    the_var = prefix+'bmod_t89_gsm'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'B!S!U'+strupcase(model)+'!R!N', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_type', 'car', $
        'coord_labels', coord_labels, $
        'model', 'T89', $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end

pro global_efield_prepare_bmod_t89_sm, probe=mission_probe, project=project, filename=cdf_file

    model = 't89'
    coord_name = 'sm'
    coord_labels = ['x','y','z']

    prefix = project[mission_probe].prefix
    var_suffix = 'bmod_'+model+'_gsm'
    global_efield_load_data, var_suffix, probe=mission_probe
    get_data, prefix+var_suffix, times, b0gsm, limits=lim
    data = cotran(b0gsm, times, 'gsm2sm')
    the_var = prefix+'bmod_t89_sm'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'B!S!U'+strupcase(model)+'!R!N', $
        'unit', 'nT', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', coord_labels, $
        'model', 'T89', $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end

pro global_efield_prepare_b_sm, probe=mission_probe, project=project, filename=cdf_file

    prefix = project[mission_probe].prefix
    var_suffix = 'b_gsm'
    global_efield_load_data, var_suffix, probe=mission_probe
    bgsm_var = prefix+var_suffix
    get_data, bgsm_var, time, bgsm, limits=lim

    data = cotran(bgsm, time, 'gsm2sm')
    the_var = prefix+'b_sm'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', 'nT', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', ['x','y','z'], $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end

pro global_efield_prepare_r_sm, probe=mission_probe, project=project, filename=cdf_file

    prefix = project[mission_probe].prefix
    var_suffix = 'r_gsm'
    global_efield_load_data, var_suffix, probe=mission_probe
    rgsm_var = prefix+var_suffix
    get_data, rgsm_var, time, rgsm, limits=lim

    data = cotran(rgsm, time, 'gsm2sm')
    the_var = prefix+'r_sm'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', ['x','y','z'], $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end


pro global_efield_prepare_u_sm, probe=mission_probe, project=project, filename=cdf_file

    prefix = project[mission_probe].prefix
    var_suffix = 'u_gsm'
    global_efield_load_data, var_suffix, probe=mission_probe
    ugsm_var = prefix+var_suffix
    get_data, ugsm_var, time, ugsm, limits=lim

    data = cotran(ugsm, time, 'gsm2sm')
    the_var = prefix+'u_sm'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'U!S!Uion!N!R', $
        'unit', 'km/s', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', ['x','y','z'], $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end

pro global_efield_prepare_e_sm, probe=mission_probe, project=project, filename=cdf_file

    prefix = project[mission_probe].prefix
    var_suffix = 'e_gsm'
    global_efield_load_data, var_suffix, probe=mission_probe
    ugsm_var = prefix+var_suffix
    get_data, ugsm_var, time, ugsm, limits=lim

    data = cotran(ugsm, time, 'gsm2sm')
    the_var = prefix+'e_sm'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'E', $
        'unit', 'mV/m', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', ['x','y','z'], $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end


pro global_efield_prepare_suxb_sm, probe=mission_probe, project=project, filename=cdf_file

    coord = 'sm'
;    b_suffix = 'b_'+coord
;    bmod_suffix = 'bmod_t89_'+coord
;    e_suffix = 'euxb_'+coord
;    global_efield_load_data, b_suffix, probe=mission_probe
;    global_efield_load_data, bmod_suffix, probe=mission_probe
;    global_efield_load_data, e_suffix, probe=mission_probe
;
;    prefix = project[mission_probe].prefix
;    e_var = prefix+e_suffix
;    b_var = prefix+b_suffix
;    bmod_var = prefix+bmod_suffix
;    sys_subtract, b_var, bmod_var, to=b_var
;    get_data, b_var, time, data
;    index = where(snorm(data) ge 5e3, count)
;    if count gt 0 then data[index,*] = !values.f_nan
;    store_data, b_var, time, data


    b_suffix = 'b_'+coord
    e_suffix = 'euxb_'+coord
    global_efield_load_data, b_suffix, probe=mission_probe
    global_efield_load_data, e_suffix, probe=mission_probe

    prefix = project[mission_probe].prefix
    e_var = prefix+e_suffix
    get_data, e_var, time, edata, limits=lim
    b_var = prefix+b_suffix
    get_data, b_var, time, bdata
    the_var = prefix+'suxb_'+coord
    data = spoynt(edata, bdata)
    cdf_save_var, the_var, value=data, filename=cdf_file
    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'S!S!UUxB!N!R', $
        'period_range', [min_period,max_period], $
        'unit', 'mW/m!U2!N', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', ['x','y','z'], $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end


; In ideal MHD, E = -UxB, where U is the plasma velocity.
pro global_efield_prepare_euxb_sm, probe=mission_probe, project=project, filename=cdf_file

    coord = 'sm'
    u_suffix = 'u_'+coord
    b_suffix = 'b_'+coord
    global_efield_load_data, u_suffix, probe=mission_probe
    global_efield_load_data, b_suffix, probe=mission_probe

    prefix = project[mission_probe].prefix
    u_var = prefix+u_suffix
    get_data, u_var, time, uvec, limits=lim
    b_var = prefix+b_suffix
    get_data, b_var, time, bvec
    ; E = -UxB, U in km/s, B in nT, E in mV/m.
    data = -vec_cross(uvec, bvec)*1e-3
    the_var = prefix+'euxb_'+coord
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'E!S!UUxB!N!R', $
        'unit', 'mV/m', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', ['x','y','z'], $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file
end

pro global_efield_prepare_pf_sm_norm, probe=mission_probe, project=project, filename=cdf_file

    pf_gsm_suffix = 'pf_gsm_norm'
    global_efield_load_data, pf_gsm_suffix, probe=mission_probe

    prefix = project[mission_probe].prefix
    pf_gsm_var = prefix+pf_gsm_suffix
    get_data, pf_gsm_var, time, pf_gsm, limits=lim
    pf_sm_suffix = 'pf_sm_norm'
    pf_sm_var = prefix+pf_sm_suffix
    data = cotran(pf_gsm, time, 'gsm2sm')
    cdf_save_var, pf_sm_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'vector', $
        'short_name', 'S', $
        'unit', 'mW/m!U2!N', $
        'coord', 'SM', $
        'coord_type', 'car', $
        'coord_labels', ['x','y','z'], $
        'depend_0', ut_name)
    cdf_save_setting, setting, varname=pf_sm_var, filename=cdf_file
end


; In MHD, pressure = ele+ion presure, P_s = n*kT, in eV/cm^3.
pro global_efield_prepare_p, probe=mission_probe, project=project, filename=cdf_file

    ele_suffix = 'ele_t'
    ion_suffix = 'ion_t'
    global_efield_load_data, ele_suffix, probe=mission_probe
    global_efield_load_data, ion_suffix, probe=mission_probe

    prefix = project[mission_probe].prefix
    get_data, prefix+ele_suffix, time, ele_temp, limits=lim
    get_data, prefix+ion_suffix, time, ion_temp
    temp = ele_temp+ion_temp
    density = 'ele_n'
    global_efield_load_data, density, probe=mission_probe
    get_data, prefix+density, time, ele_n
    index = where(abs(ele_n) ge 1e3, count)
    if count eq 0 then ele_n[index] = !values.f_nan
    data = ele_n*temp
    the_var = prefix+'p'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'scalar', $
        'short_name', 'P', $
        'unit', 'eV/cm!U3!N', $
        'depend_0', ut_name, $
        'ylog', 1)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end


; B^2/2mu0, in eV/cm^3.
pro global_efield_prepare_pmag, probe=mission_probe, project=project, filename=cdf_file

    b_suffix = 'b_gsm'
    global_efield_load_data, b_suffix, probe=mission_probe
    prefix = project[mission_probe].prefix
    get_data, prefix+b_suffix, time, bgsm, limits=lim

    mu0 = 4*!dpi*1e-7
    coef1 = 1e-18/(2*mu0)   ; 1 nT, to SI.
    coef2 = 1.6e-19*1e6     ; eV/cm^2 to SI.
    data = total(bgsm^2,2)*(coef1/coef2)

    the_var = prefix+'pmag'
    cdf_save_var, the_var, value=data, filename=cdf_file

    ut_name = lim.depend_0
    setting = dictionary($
        'display_type', 'scalar', $
        'short_name', 'P!UB!N', $
        'unit', 'eV/cm!U3!N', $
        'depend_0', ut_name, $
        'ylog', 1)
    cdf_save_setting, setting, varname=the_var, filename=cdf_file

end


pro global_efield_load_data, data_name, probe=mission_probe, project=project, _extra=ex

    if n_elements(data_name) eq 0 then message, 'No data_name ...'
    if n_elements(mission_probe) eq 0 then message, 'No mission_probe'
    if n_elements(project) eq 0 then project = global_efield_load_project()
    index = where(project.all_mission_probes eq mission_probe, count)
    if count eq 0 then message, 'Invalid mission_probe: '+mission_probe+' ...'

    file_key = 'combined_file_suffix'
    if ~project[mission_probe].haskey(file_key) then $
        global_efield_prepare_primitive_data, 'orbit', probe=mission_probe, project=project
    cdf_file = join_path([project.data_dir,(project[mission_probe])[file_key]])
    if file_test(cdf_file) eq 0 then begin
        lprmsg, 'Combined_file does not exist, resttting ...'
        project[mission_probe].remove, file_key
        update_project, project
        global_efield_load_data, data_name, probe=mission_probe, project=project, _extra=ex
    endif


    prefix = project[mission_probe].prefix
    var = prefix+data_name
    if ~cdf_has_var(var, filename=cdf_file) then begin
    ;---Need to load data to CDF first.
        case data_name of
            'b_gsm': routine_suffix = 'bfield'
            'e_gsm': routine_suffix = 'efield'
            'u_gsm': routine_suffix = 'ion_vel'
            'ele_n': routine_suffix = 'density'
            'r_gsm': routine_suffix = 'orbit'
            'ele_t': routine_suffix = 'ele_temp'
            'ion_t': routine_suffix = 'ion_temp'
            'pf_gsm_norm': routine_suffix = 'pflux'
            else: routine_suffix = !null
        endcase

        if n_elements(routine_suffix) ne 0 then begin
        ;---Load primitive quantities.
            global_efield_prepare_primitive_data, routine_suffix, probe=mission_probe, project=project
        endif else begin
        ;---Load derived quantities.
            routine_name = project.name+'_prepare_'+data_name
            call_procedure, routine_name, probe=mission_probe, project=project, filename=cdf_file, _extra=ex
        endelse
    endif

    cdf_load_var, var, filename=cdf_file

end

mission_probes = ['rbsp'+letters('b'),'th'+letters('e')]
var_suffixes = ['r_sm','b_sm','bmod_t89_sm']
mission_probes = ['th'+letters('e'),'polar','rbsp'+letters('b')]
var_suffixes = ['e_gsm']

;mission_probes = 'c'+['1','2','3','4']
;var_suffixes = ['b_sm','bmod_t89_sm']
;mission_probes = 'th'+letters('e')
;var_suffixes = ['e_gsm']

mission_probes = ['polar','rbsp'+letters('b'),'th'+letters('e'),'c'+['1','2','3','4']]
var_suffixes = ['bmod_t89_gsm','b_gsm']

mission_probes = 'rbspb'
var_suffixes = 'pf_sm_norm'

foreach mission_probe, mission_probes do foreach var_suffix,var_suffixes do global_efield_load_data, var_suffix, probe=mission_probe
end