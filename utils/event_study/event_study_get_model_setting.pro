
function event_study_get_model_setting, event_info

    key = 'model_setting'
    if ~event_info.haskey(key) then begin
        model_setting = dictionary($
            'external_model', 't89', $
            'models', ['t89','t96','t01','t04s'], $
            't89_par', 2, $
            'igrf', 0, $
            'refine', 1, $
            'direction', 'north' )
        internal_model = (model_setting['igrf'] eq 0)? 'dipole': 'igrf'
        model_setting['internal_model'] = internal_model
        event_info[key] = model_setting
    endif
    return, event_info[key]

end