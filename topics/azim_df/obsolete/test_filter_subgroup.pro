;+
; Test to filter subgroups using azim_df_analyze_subgroup.
; 
; Looks good. Run one time and done.
;-

project = azim_df_load_project()
search_names = ['within_15Re','beyond_15Re']
routines = 'azim_df_search_'+['post_midn','pre_midn','around_midn']+'_events'

good_list = join_path([project.data_dir,'azim_df_analyze_subgroup_good_ones.txt'])
bad_list = join_path([project.data_dir,'azim_df_analyze_subgroup_bad_ones.txt'])
;foreach file, [good_list,bad_list] do begin
;    if file_test(file) eq 1 then file_delete, file
;    ftouch, file
;endforeach

foreach search_name, search_names do begin
    foreach routine, routines do begin
        candidates = call_function(routine, search_step='subgroup')
        foreach candidate, candidates do begin
            if candidate.search_type ne search_name then continue
            stop    ; put a stop here to prevent unintensional runs.
            subgroup = azim_df_analyze_subgroup(candidate, project=project)
            file = (n_elements(subgroup) eq 0)? bad_list: good_list
            azim_df_subgroup_write_file, candidate, filename=file
        endforeach
    endforeach
endforeach

end