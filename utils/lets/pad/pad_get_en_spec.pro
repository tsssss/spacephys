;+
; Get en spec from pad_var.
;-
function pad_get_en_spec, pad_var=pad_var, pitch_angle_range=pitch_angle_range, $
    var_info=var_info

    if n_elements(var_info) eq 0 then var_info = streplace(pad_var,'pad','en_spec')
    pad_fluxs = get_var_data(pad_var, times=times, settings=settings)

    pa_centers = settings.pa_centers
    en_centers = settings.en_centers

    if n_elements(pitch_angle_range) ne 2 then pitch_angle_range = [0d,180]
    pa_index = where_pro(pa_centers, '[]', pitch_angle_range, count=count)
    if count eq 0 then begin
        errmsg = 'Invalid pitch angle range ...'
        return, retval
    endif

    en_fluxs = total(pad_fluxs[*,pa_index,*],2,nan=1)/count
    store_data, var_info, times, en_fluxs, en_centers
    unit = settings.unit
    species = settings.species
    species_str = species
    energy_unit = settings.en_unit
    time_range = settings.requested_time_range
    add_setting, var_info, smart=1, {$
        requested_time_range: time_range, $
        display_type: 'spec', $
        unit: unit, $
        species: species, $
        species_name: species_str, $
        ytitle: 'Energy!C('+energy_unit+')', $
        subytitle: species_str, $
        ylog: 1, $
        zlog: 1, $
        short_name: ''}

    return, var_info

end