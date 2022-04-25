;+
; Return the nominal energy bins for a given probe and species.
;-

function ml_rbsp_hope_read_energy_bins, probe=probe, species=species

    all_species = rbsp_hope_species()
    index = where(species eq all_species, count)
    if count eq 0 then begin
        errmsg = 'Invalid species '+species+' ...'
        return, !null
    endif
    charge_type = (species[0] eq 'e')? 'ele': 'ion'

    energy_bin_file = join_path([srootdir(),'rbsp_hope_energy_bins.cdf'])
    if file_test(energy_bin_file) eq 0 then begin
        cdf_touch, energy_bin_file
    endif

    prefix = 'rbsp'+probe+'_'
    var = prefix+'energy_'+charge_type
    if ~cdf_has_var(var, filename=energy_bin_file) then begin
        start_day = time_double('2013-01-01')
        secofday = constant('secofday')
        go = 1
        energy_bin_var = (charge_type eq 'ele')? 'HOPE_ENERGY_Ele': 'HOPE_ENERGY_Ion'
        lowest_energy_bin = (charge_type eq 'ele')? 15d: 1d
        settings = dictionary($
            'charge_type', charge_type, $
            'unit', 'eV' )
        while go do begin
            time_range = start_day+[0,secofday]
            files = rbsp_load_hope(time_range, probe=probe, id='l2%sector')
            if file_test(files[0]) eq 0 then begin
                start_day += secofday
                continue
            endif
            energy_bins = cdf_read_var(energy_bin_var, filename=files[0])
            lowest_energies = energy_bins[*,0]
            index = where(lowest_energies le lowest_energy_bin, count)
            if count ne 0 then begin
                energy_bins = reform(energy_bins[index[0],*])
                cdf_save_var, var, value=energy_bins, filename=energy_bin_file
                cdf_save_setting, settings, varname=var, filename=energy_bin_file
                break
            endif
            start_day += secofday
        endwhile
    endif


    return, cdf_read_var(var, filename=energy_bin_file)

end

tmp = ml_rbsp_hope_read_energy_bins(probe='b', species='e')
end