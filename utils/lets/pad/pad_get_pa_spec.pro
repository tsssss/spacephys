;+
; Get pa spec from pad_var.
;-
function pad_get_pa_spec, pad_var=pad_var, energy_range=energy_range, var_info=var_info

    if n_elements(var_info) eq 0 then var_info = streplace(pad_var,'pad','pa_spec')
    pad_fluxs = get_var_data(pad_var, times=times, settings=settings)

    pa_centers = settings.pa_centers
    en_centers = settings.en_centers

    if n_elements(energy_range) ne 2 then energy_range = minmax(en_centers)
    en_index = where_pro(en_centers, '[]', energy_range, count=count)
    if count eq 0 then begin
        errmsg = 'Invalid energy range ...'
        return, retval
    endif

    pa_fluxs = total(pad_fluxs[*,*,en_index],3,nan=1)/count
    store_data, var_info, times, pa_fluxs, pa_centers
    unit = settings.unit
    species = settings.species
    species_str = species
    pa_unit = settings.pa_unit
    time_range = settings.requested_time_range
    add_setting, var_info, smart=1, {$
        requested_time_range: time_range, $
        display_type: 'spec', $
        unit: unit, $
        species: species, $
        species_name: species_str, $
        ytitle: 'PA!C('+pa_unit+')', $
        subytitle: species_str, $
        ylog: 0, $
        zlog: 1, $
        yrange: [0,180], $
        ytickv: [0,90,180], $
        yticks: 2, $
        yminor: 3, $
        short_name: ''}

    return, var_info

end