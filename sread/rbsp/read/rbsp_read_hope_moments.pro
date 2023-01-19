;+
; Read RBSP HOPE moments.
;-

function rbsp_read_hope_moments, input_time_range, probe=probe, errmsg=errmsg, $
    species=species, coord=coord
    

    errmsg = ''
    retval = ''
    prefix = 'rbsp'+probe+'_'

    if n_elements(species) eq 0 then species = 'p'
    index = where(species eq rbsp_hope_species(), count)
    if count eq 0 then begin
        errmsg = 'Invalid species: '+species+' ...'
        return, retval
    endif
    species_name = rbsp_hope_species_name(species)
    if n_elements(coord) eq 0 then coord = 'gsm'

    vars = prefix+species+'_'+['n','t',['vbulk','nflux','eflux','enthalpy']+'_'+coord]
    
    time_range = time_double(input_time_range)
    files = rbsp_load_hope_moments(time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval

    var_list = list()

    default_coord = 'gse'
    vec_default_vars = prefix+species+'_'+['vbulk','nflux','eflux','enthalpy']+'_'+default_coord
    scalar_vars = prefix+species+'_'+['density','t_avg']
    in_vars = [scalar_vars,vec_default_vars]
    out_vars = [vars[0:1],vec_default_vars]
    time_var = (species eq 'e')? 'ut_ele': 'ut_ion'
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'out_vars', out_vars, $
        'time_var_name', time_var, $
        'time_var_type', 'unix' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    ; Conver to the wanted coord.
    if coord ne default_coord then begin
        vec_coord_vars = prefix+species+'_'+['vbulk','nflux','eflux','enthalpy']+'_'+coord
        short_names = ['V','F',tex2str('Gamma'),'H']
        foreach vec_default_var, vec_default_vars, var_id do begin
            vec_coord_var = vec_coord_vars[var_id]
            get_data, vec_default_var, times, vec_default, limits=lim
            vec_coord = cotran(vec_default, times, default_coord+'2'+coord, probe=probe)
            store_data, vec_coord_var, times, vec_coord, limits=lim

            short_name = species_name+' '+short_names[var_id]
            vatt = cdf_read_setting(vec_default_var, filename=files[0])
            unit = vatt.units
            add_setting, vec_coord_var, /smart, {$
                display_type: 'vector', $
                unit: unit, $
                short_name: short_name, $
                coord: strupcase(coord), $
                coord_labels: ['x','y','z'], $
                colors: constant('rgb') }
        endforeach
    endif

    short_names = ['N','T']
    foreach scalar_var, scalar_vars, var_id do begin
        short_name = species_name+' '+short_names[var_id]
        vatt = cdf_read_setting(vec_default_var, filename=files[0])
        unit = vatt.units

        add_setting, scalar_var, /smart, {$
            display_type: 'scalar', $
            unit: unit, $
            short_name: short_name, $
            ylog: 1 }
    endforeach

    return, vars

end


time_range = ['2013-05-01','2013-05-02']
probe = 'b'
species = 'p'
vars = rbsp_read_hope_moments(time_range, probe=probe, species=species)
end