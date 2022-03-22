;+
; Read HOPE moments.
;-

pro rbsp_load_hope_moments_gen_file, input_time_range, probe=probe, $
    errmsg=errmsg, filename=cdf_file


    secofday = constant('secofday')
    date = time_double(input_time_range[0])
    date = date-(date mod secofday)
    time_range = date+[0,secofday]
    ion_energy_range = [30,1e6]
    electron_energy_range = [200,1e6]
    species = ['e','p','o','he']
    coord = 'gse'

    rbsp_calc_hope_moments, time_range, probe=probe, errmsg=errmsg, $
        species=species, coord=coord, $
        ion_energy_range=ion_energy_range, electron_energy_range=electron_energy_range

    if errmsg ne '' then return

;---Save data to cdf file.
    if n_elements(cdf_file) eq 0 then begin
        errmsg = handle_error('No output CDF file ...')
        return
    endif
    if file_test(cdf_file) ne 0 then file_delete, cdf_file
    gatt = dictionary($
        'title', 'RBSP HOPE moments, calculated based on l2 data', $
        'text', 'Calculated by Sheng Tian at the University of Minnesota, email:tianx138@umn.edu' )
    cdf_save_setting, gatt, filename=cdf_file

    prefix = 'rbsp'+probe+'_'
    vatt = dictionary($
        'FIELDNAM', 'Unix time', $
        'UNITS', 'sec', $
        'VAR_TYPE', 'support_data' )
    foreach the_species, species do begin
        charge_type = (the_species eq 'e')? 'ele': 'ion'
        time_var = 'ut_'+charge_type
        if ~cdf_has_var(time_var, filename=cdf_file) then begin
            moms = get_var_data(prefix+the_species+'_moms')
            times = moms.ut
            cdf_save_var, time_var, value=times, filename=cdf_file, cdf_type='CDF_DOUBLE'
            cdf_save_setting, vatt, varname=time_var, filename=cdf_file
        endif

        vars = prefix+the_species+'_'+['density','t_avg',$
            ['vbulk','nflux','eflux','enthalpy']+'_'+coord]
        foreach var, vars do begin
            data = float(get_var_data(var, limits=lim))
            ndim = size(data, n_dimensions=1)
;            if ndim gt 1 then data = transpose(data)
            settings = dictionary()
            settings['DEPEND_0'] = time_var
            settings['VAR_TYPE'] = 'data'
            settings['UNITS'] = lim.unit
            settings['DISPLAY_TYPE'] = lim.display_type
            cdf_save_var, var, value=data, filename=cdf_file
            cdf_save_setting, settings, varname=var, filename=cdf_file
        endforeach
    endforeach

end

time_range = ['2013-06-07','2013-06-08']
probe = 'a'
file = join_path([homedir(),'test_hope_moments.cdf'])
file_delete, file, allow_nonexistent=1
rbsp_load_hope_moments_gen_file, time_range, probe=probe, filename=file
end
