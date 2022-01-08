;+
; A wrapper to load events azimuthally away from midnight.
;-

function azim_df_load_azim_events, project=project

    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_load_coherent_events(project=project)
    index_info = azim_df_subgroup_analyze_direction(project=project)
    index = index_info.azim_away.index
    if n_elements(index) eq 0 then return, !null
    return, events[index]

end
