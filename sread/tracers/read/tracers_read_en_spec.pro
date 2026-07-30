

function tracers_read_en_spec, input_time_range, probe=probe, errmsg=errmsg, get_name=get_name, suffix=suffix, species=species, _extra=ex

    errmsg = ''

    sources = ['iowa','cdaweb']
    species_info = dictionary( $
        'e', 'ele', $
        'i', 'ion' )
    the_species = strlowcase(species)
    if not species_info.haskey(the_species) then begin
        errmsg = 'Invalid species: '+species+' ...'
        return, !null
    endif
    foreach source, sources do begin
        func_name = 'tracers_read_'+species_info[the_species]+'_en_spec_'+source
        retval = call_function(func_name, input_time_range, probe=probe, errmsg=errmsg, get_name=get_name, suffix='', _extra=ex)
        if errmsg eq '' then return, retval
    endforeach


end