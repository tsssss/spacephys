;+
; Read s/c position.
;
; time.
; probe=.
; coord=.
;-

pro azim_dp_read_orbit, time, probe=probe, errmsg=errmsg, coord=coord

    prefix = probe+'_'
    if n_elements(coord) eq 0 then coord = 'sm'
    azim_dp_read_level1_data, time, probe=probe, datatype='r_sm', errmsg=errmsg
    if errmsg ne '' then return

    r_sm_var = prefix+'r_sm'
    r_sm = get_var_data(r_sm_var, times=times)
    index = where(finite(r_sm[*,0]), count)
    if count eq 0 then begin
        errmsg = 'No data ...'
        del_data, r_sm_var
        return
    endif
    
    the_coord = strlowcase(coord)
    r_coord_var = prefix+'r_'+the_coord
    if the_coord ne 'sm' then begin
        r_coord = cotran(r_sm, times, 'sm2'+the_coord)
        store_data, r_coord_var, times, r_coord
    endif
    add_setting, r_coord_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'Re', $
        'short_name', 'R', $
        'coord', strupcase(the_coord), $
        'coord_labels', constant('xyz') )

end

time = time_double(['2019-08-01','2019-09-01'])
azim_dp_read_orbit, time, probe='tha', coord='sm'
end
