;+
; Read solar wind B in given coord.
;-

function ml_omni_read_sw_b, input_time_range, get_name=get_name, coord=input_coord

    errmsg = ''
    retval = ''

;---Check input coord.
    if n_elements(input_coord) eq 0 then input_coord = 'gse'
    coord = strlowcase(input_coord[0])

;---Return if only name is needed.
    prefix = 'omni_'
    coord_var = prefix+'sw_b_'+coord
    if keyword_set(get_name) then return, coord_var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)

    gse_var = ml_omni_read_param_var(input_time_range, var='sw_b_gse', errmsg=errmsg)
    if errmsg ne '' then return, retval

    if coord ne 'gse' then begin
        get_data, gse_var, times, vec_gse
        vec_coord = cotran(vec_gse, times, 'gse2'+coord)
        store_data, coord_var, times, vec_coord
    endif
    add_setting, coord_var, smart=1, dictionary($
        'display_type', 'vector', $
        'unit', 'nT', $
        'short_name', 'SW B', $
        'coord', strupcase(coord), $
        'coord_labels', constant('xyz') )
    return, coord_var

end


time_range = ['2008-01-17','2008-01-18']
var = ml_omni_read_sw_b(time_range, coord='sm')
end