;+
; Read R in Re.
; 
; resolution=. A number for cadence, 5*60 or 60 sec.
;-

function ml_rbsp_read_pos, input_time_range, probe=probe, $
    coord=input_coord, errmsg=errmsg, get_name=get_name, resolution=resolution

    errmsg = ''
    retval = ''

;---Check coord.
    if n_elements(input_coord) eq 0 then input_coord = 'sm'
    coord = strlowcase(input_coord[0])

;---Check input probe.
    probes = rbsp_probes()
    index = where(probes eq probe, count)
    if count eq 0 then begin
        errmsg = 'probe is unknown ...'
        return, retval
    endif
    prefix = 'rbsp'+probe+'_'

;---Return if only name is needed.
    r_coord_var = prefix+'r_'+coord
    if keyword_set(get_name) then return, r_coord_var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)

;---Read files.
    files = ml_rbsp_load_orbit_var(time_range, probe=probe, errmsg=errmsg, resolution=60)
    if errmsg ne '' then return, retval

;---Read vars.
    var_list = list()
    r_sm_var = prefix+'r_sm'
    var_list.add, dictionary($
        'in_vars', r_sm_var, $
        'time_var_name', 'unix_time', $
        'time_var_type', 'unix' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    if coord ne 'sm' then begin
        get_data, r_sm_var, times, r_sm
        r_coord = cotran(r_sm, times, 'sm2'+coord)
        store_data, r_coord_var, times, r_coord
    endif
    add_setting, r_coord_var, smart=1, dictionary($
        'requested_time_range', time_range, $
        'display_type', 'vector', $
        'unit', 'Re', $
        'short_name', 'R', $
        'coord', strupcase(coord), $
        'coord_labels', constant('xyz') )
    return, r_coord_var
end

time_range = ['2013-06-07','2013-06-08']
probe = 'a'
var = ml_rbsp_read_pos(time_range, probe=probe)
end