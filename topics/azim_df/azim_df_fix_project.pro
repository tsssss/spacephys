;+
; The tplot file saving project can be damaged.
; This program fix the tplot file using existing data files.
;-

pro azim_df_fix_project

;---Fix the basic project.
    project_name = 'azim_df'
    project = init_project(project_name)
    project = azim_df_load_project()
    
;---Fix the settings and results of candidate_search.
    azim_df_search_candidate, project=project, /fix

end