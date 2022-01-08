; Init project.
project = azim_df_load_project()

; Generate diagnostic plot for primitive data (r_gsm and b_gsm).
azim_df_gen_diagnostic_plot, 'primitive_data', project=project

; Search candidate based on orbit data, omni data, and avaiability of b field data.
candidates = azim_df_search_candidate(project=project)
azim_df_gen_diagnostic_plot, 'search_candidate_search_roi', project=project
azim_df_gen_diagnostic_plot, 'search_candidate_search_triad', project=project
azim_df_gen_diagnostic_plot, 'search_candidate', project=project

;; This is for searching large DF events.
;region_names = ['pre','post']+'_midn'
;search_names = ['within_15Re','beyond_15Re']
;foreach search_name, search_names do foreach region_name, region_names do $
;    azim_df_gen_diagnostic_plot, 'search_event_search_large_df', project=project, $
;    region_name=region_name, search_name=search_name
;; This is for searching DF groups.
;azim_df_gen_diagnostic_plot, 'search_event_search_df_group', project=project

; Then I manually modified the results of event search to:
;   1. Add small DFs that are missing in the search.
;   2. Correct for DFs that are mis-grouped by the program.
;   All these corrections are using the existing DFs in the df_search.








; Find candidates in 2012-2017.
if ~project.haskey('done_search_candidate') then azim_df_search_candidate, project=project
; Find themis candidates in 2007-2009.
if ~project.haskey('done_search_candidate_themis') then azim_df_search_candidate_themis, project=project

; Load data for selected storms.
candidate = project.candidate
foreach event_id, candidate.keys() do begin
    tinfo = candidate[event_id]
    the_key = 'data_file_suffix'
    if ~tinfo.haskey(the_key) then tinfo[the_key] = event_id+'_all_data.tplot'
    data_file = join_path([project.data_dir,tinfo[the_key]])
    if file_test(data_file) eq 0 then begin
        event_time_range = tinfo.time_range
        azim_df_load_basic_data, event_time_range, event_id=event_id, project=project, reset=reset
    endif
endforeach
update_project, project
stop



; Need to manually decide the time ranges for specific event.
; Save them to <data_dir>/event_list.txt.

; Use azim_df_load_basic_data to load r_sm,mlt,mlat,db_tilt,theta, and
; pre/post_midn_triad_flag, pre/post_midn_data_flag.

; Use thata to calculate the time_lag, ref_time, triad_probes.

end
