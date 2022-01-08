;+
; A wrapper to load the globally coherent events.
;-

function azim_df_load_coherent_events, project=project
    events = azim_df_subgroup_analyze_overlap(project=project)
    events = azim_df_filter_triads(events, project=project)
    events = azim_df_filter_v2d(events, project=project)
    return, events
end
