;+
; Read B field.
;
; time.
; probe=.
; coord=.
;-

pro azim_dp_read_bfield, time, probe=probe, errmsg=errmsg, coord=coord

    prefix = probe+'_'
    if n_elements(coord) eq 0 then coord = 'sm'
    azim_dp_read_level1_data, time, probe=probe, datatype='b_sm', errmsg=errmsg
    if errmsg ne '' then return

    b_sm_var = prefix+'b_sm'
    b_sm = get_var_data(b_sm_var, times=times)
    index = where(finite(b_sm[*,0]), count)
    if count eq 0 then begin
        errmsg = 'No data ...'
        del_data, b_sm_var
        return
    endif

    the_coord = strlowcase(coord)
    b_coord_var = prefix+'b_'+the_coord
    if the_coord ne 'sm' then begin
        b_coord = cotran(b_sm, times, 'sm2'+the_coord)
        store_data, b_coord_var, times, b_coord
    endif
    add_setting, b_coord_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'nT', $
        'short_name', 'B', $
        'coord', strupcase(the_coord), $
        'coord_labels', constant('xyz') )

end

time = time_double(['2019-08-01','2019-09-01'])
azim_dp_read_bfield, time, probe='tha', coord='sm'
end
