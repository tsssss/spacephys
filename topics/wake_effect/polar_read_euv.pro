;+
; Read Euv for Polar.
;-

pro polar_read_euv, time, probe=probe

    spin_period = polar_info('spin_period')
    data_rate = polar_info('v_uvw_data_rate')

    polar_read_ebv, time, id='v_uvw', errmsg=errmsg
    if errmsg ne '' then return
    v_var = 'po_v_uvw'
    uniform_time, v_var, data_rate

;---Convert V to mV/m, use the nominal boom length, ignore shortening factors.
    get_data, v_var, times, euv
    boom_lengths = polar_info('boom_length')
    euv = euv[*,0:1]
    for ii=0,1 do euv[*,ii] = euv[*,ii]/boom_lengths[ii]*1e3    ; V to mV/m.

    nspin = 2d  ; width for calculate the DC offset.
    width = spin_period*nspin/data_rate
    for ii=0, 1 do begin
        dc_offset = smooth(euv[*,ii], width, /edge_truncate, /nan)
        dc_offset = smooth(dc_offset, width, /edge_truncate, /nan)
        euv[*,ii] -= dc_offset
    endfor
    e_var = 'po_e_uv'
    store_data, e_var, times, euv
    add_setting, e_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E', $
        coord: 'UVW', $
        coord_labels: ['u','v'], $
        colors: sgcolor(['red','blue'])}

end

time = time_double('1998-09-25')+[0,86400d]
polar_read_euv, time
end
