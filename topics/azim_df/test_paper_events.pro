;+
; Do all needed analysis for selected events.
;-

project = azim_df_load_project()
file = join_path([project.data_dir,'azim_df_paper_events.txt'])

events = azim_df_subgroup_read_file(file)
foreach candidate, events do begin
    event = azim_df_filter_coherent_df(candidate, project=project)
;    if n_elements(event) eq 0 then stop
;    azim_df_event_survey_plot2, event, project=project, dirname='paper_events'
endforeach

end