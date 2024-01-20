function inherit_setting, settings, id=id

    if n_elements(id) eq 0 then id = 'basic'

    if id eq 'basic' then keys = ['requested_time_range',$
        'mission_probe']
    
    new_settings = dictionary()
    foreach key, keys do begin
        if ~settings.haskey(key) then continue
        new_settings[key] = settings[key]
    endforeach

    return, new_settings

end