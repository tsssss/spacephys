;+
; Read HOPE flux for a given species, energy, and probe.
;-

function ml_rbsp_hope_read_flux, input_time_range, probe=probe, $
    species=input_species, energy=input_energy, to=out_var, errmsg=errmsg, get_name=get_name, resolution=resolution

    errmsg = ''
    retval = ''

;---Check input energy.
    if n_elements(input_energy) ne 1 then begin
        errmsg = 'Invalid input_energy ...'
        return, retval
    endif
    energy = float(input_energy)

;---Check input species and probe.
    en_spec_var = ml_rbsp_hope_read_en_spec(probe=probe, species=input_species, get_name=1, resolution=resolution)
    if en_spec_var eq '' then begin
        errmsg = 'Invalid input probe, or species ...'
        return, retval
    endif
    prefix = 'rbsp'+probe+'_'

;---Get the actual energy.
    species = strlowcase(input_species[0])
    energy_bins = ml_rbsp_hope_read_energy_bins(probe=probe, species=species)
    energy_range = minmax(energy_bins)
    if product(energy_range-energy) gt 0 then begin
        errmsg = 'input_energy is invalid, should be in '+strjoin(string(round(energy_range),format='(I0)'),',')+' eV'
        return, retval
    endif
    tmp = min(energy_bins-energy, abs=1, energy_index)
    energy = energy_bins[energy_index]
    energy_str = string(round(energy),format='(I0)')

;---Return if only name is needed.
    flux_var = en_spec_var+'_'+energy_str+'eV'
    if keyword_set(get_name) then return, flux_var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)

    time_step = ml_time_step()
    if check_if_update(en_spec_var, time_range, dtime=time_step) then begin
        en_spec_var = ml_rbsp_hope_read_en_spec(time_range, probe=probe, species=input_species, resolution=resolution)
    endif

    get_data, en_spec_var, times, en_specs
    flux = en_specs[*,energy_index]
    store_data, flux_var, times, flux
    add_setting, flux_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'unit', 'eV', $
        'short_name', 'F', $
        'ylog', 1, $
        'yrange', [1,1e9], $
        'species', species )
    return, flux_var

end

time_range = ['2013-04-30','2013-05-02']
probe = 'b'
species = 'p'
energy = 100
flux_var = ml_rbsp_hope_read_flux(time_range, probe=probe, species=species, energy=energy)
end
