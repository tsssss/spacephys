;+
; Read Euv for RBSP.
;-

pro rbsp_read_euv, time, probe=probe, errmsg=errmsg

    errmsg = ''
    spin_period = rbsp_info('spin_period')
    data_rate = rbsp_info('v_uvw_data_rate')

    pre0 = 'rbsp'+probe+'_'

    rbsp_read_efield, time, probe=probe, resolution='hires'
    e0_var = pre0+'e0_gsm'
    get_data, e0_var, times, egsm
    store_data, e0_var, times, egsm

    rbsp_read_orbit, time, probe=probe, errmsg=errmsg
    if errmsg ne '' then return
    rbsp_read_bfield, time, probe=probe, resolution='4sec', errmsg=errmsg
    if errmsg ne '' then return
    rbsp_read_quaternion, time, probe=probe, errmsg=errmsg
    if errmsg ne '' then return
    rbsp_calc_emodel, time, probe=probe, errmsg=errmsg
    if errmsg ne '' then return

    emod_var = pre0+'emod_gsm'
    get_data, emod_var, uts, e0gsm
    e0gsm = sinterpol(e0gsm, uts, times, /nan)
    e0mag = snorm(e0gsm)
    emag = snorm(egsm)
    index = where(finite(emag))
    fit_res = linfit(e0mag[index], emag[index])
    fs_coef = fit_res[1]
    degsm = egsm/fs_coef-e0gsm
    de_var = pre0+'de_gsm'
    store_data, de_var, times, degsm, limits=lims
    add_setting, de_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E', $
        coord: 'GSM', $
        coord_labels: ['x','y','z'], $
        colors: sgcolor(['red','green','blue'])}

    e_var = pre0+'e_uv'
    rbsp_gsm2uvw, de_var, e_var
    get_data, e_var, times, euv
    euv = euv[*,0:1]

    nspin = 2d  ; width for calculate the DC offset.
    width = spin_period*nspin/data_rate
    for ii=0, 1 do begin
        dc_offset = smooth(euv[*,ii], width, /edge_truncate, /nan)
        dc_offset = smooth(dc_offset, width, /edge_truncate, /nan)
        euv[*,ii] -= dc_offset
    endfor
    store_data, e_var, times, euv
    add_setting, e_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E', $
        coord: 'UVW', $
        coord_labels: ['u','v'], $
        colors: sgcolor(['red','blue'])}

end

time = time_double('2013-09-25')+[0,86400d]
rbsp_read_euv, time, probe='a'
end
