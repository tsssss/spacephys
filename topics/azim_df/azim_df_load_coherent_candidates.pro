;+
; A wrapper to load the globally coherent events.
;-

function azim_df_load_coherent_candidates, project=project
    return, azim_df_subgroup_analyze_overlap(project=project)
end
