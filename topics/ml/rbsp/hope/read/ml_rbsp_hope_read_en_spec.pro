;+
; Read HOPE en spec.
;-

function ml_rbsp_hope_read_en_spec, input_time_range, probe=probe, $
    species=input_species, to=out_var, errmsg=errmsg, get_name=get_name, resolution=resolution

    errmsg = ''
    retval = ''

;---Check input probe.
    probes = rbsp_probes()
    index = where(probes eq probe, count)
    if count eq 0 then begin
        errmsg = 'probe is unknown ...'
        return, retval
    endif
    prefix = 'rbsp'+probe+'_'

;---Check input species.
    if n_elements(input_species) eq 0 then begin
        errmsg = 'No input species ...'
        return, retval
    endif
    species = strlowcase(input_species[0])
    all_species = rbsp_hope_species()
    index = where(all_species eq species, count)
    if count eq 0 then begin
        errmsg = 'input_species is unknown, all known species are '+strjoin(all_species,',')+' ...'
        return, retval
    endif

;---Return if only name is needed.
    en_spec_var = prefix+species+'_en_spec'
    if keyword_set(get_name) then return, en_spec_var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)


;---Read files.
    files = ml_rbsp_hope_load_en_spec(time_range, probe=probe, errmsg=errmsg, resolution=resolution)
    if errmsg ne '' then return, retval

;---Read vars.
    var_list = list()
    flux_var = prefix+'hope_flux_'+species
    var_list.add, dictionary('in_vars', flux_var)
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    if n_elements(out_var) ne 0 then en_spec_var = out_var
    get_data, flux_var, times, fluxes
    energy_bins = ml_rbsp_hope_read_energy_bins(probe=probe, species=species)
    store_data, en_spec_var, times, fluxes, energy_bins
    short_name = rbsp_hope_species_name(species)
    add_setting, en_spec_var, smart=1, dictionary($
        'display_type', 'spec', $
        'unit', '#/cm!U-2!N-s-sr-keV', $
        'short_name', short_name, $
        'ylog', 1, $
        'yrange', [4,4e4], $
        'zlog', 1, $
        'species', species, $
        'requested_time_range', time_range )
    return, en_spec_var
    
end