
function azim_df_search_post_midn_events, search_step=search_step, project=project, $
    reset=reset, test_time=test_time

    project = azim_df_load_project()

;---Some common settings.
    region_name = 'post_midn'
    mlt_limit = 9.  ; hour.
    mlt_range = [0,1]*mlt_limit
    search_settings = azim_df_search_event_settings()
    foreach search_setting, search_settings do begin
        search_setting.search_roi.mlt_range = mlt_range
        search_name = search_setting.name
        foreach key, search_setting.keys() do begin
            the_setting = search_setting[key]
            if ~isa(the_setting,'dictionary') then continue
            if ~the_setting.haskey('file_suffix') then continue
            the_prefix = 'azim_df_search_event_'+region_name+'_'+search_name+'_'
            the_setting['file_suffix'] = the_prefix+key+'.txt'
        endforeach
    endforeach


    ; The full routine list.
    ;routines = 'azim_df_search_'+['roi','triad','large_df','subgroup','df_group','uniq_subgroup']
    if n_elements(search_step) eq 0 then search_step = 'uniq_subgroup'
    routine = 'azim_df_search_'+search_step
    candidates = list()
    foreach search_setting, search_settings do begin
        candidates.add, call_function(routine, search_setting, project=project, reset=reset, test_time=test_time), /extract
    endforeach

    return, candidates
end

test_time = time_double('2007-11-20/17:10')
;test_time = time_double('2014-08-28/10:20')
;test_time = time_double('2016-11-02/00:55')
;test_time = time_double('2017-03-28/01:42')
test_time = !null

reset = 1
search_step = 'uniq_subgroup'
candidates = list()
candidates.add, azim_df_search_post_midn_events(search_step=search_step, reset=reset), /extract
candidates.add, azim_df_search_pre_midn_events(search_step=search_step, reset=reset), /extract
stop


;reset = 0
;search_step = 'subgroup'
;candidates = list()
;candidates.add, azim_df_search_post_midn_events(search_step=search_step, reset=reset), /extract
;candidates.add, azim_df_search_pre_midn_events(search_step=search_step, reset=reset), /extract
;dirname = 'all_subgroups'
;azim_df_subgroup_gen_diagnostic_plot, candidates, dirname=dirname, project=project
;stop



;reset = 1
;search_steps = ['large_df','subgroup','df_group']
;foreach search_step, search_steps do begin
;    candidates = list()
;    candidates.add, azim_df_search_post_midn_events(search_step=search_step, reset=reset), /extract
;    candidates.add, azim_df_search_pre_midn_events(search_step=search_step, reset=reset), /extract
;endforeach
;stop

events = list()
selected_index = list()
rejected_index = list()
log_file = -1
lprmsg, 'Start from '+string(candidates.length,format='(I0)')+' ...', log_file
foreach candidate, candidates, candidate_id do begin
    lprmsg, 'Processing '+string(candidate_id,format='(I0)')+' ...', log_file
    event = azim_df_filter_df_group(candidate, project=project, log_file=log_file)
    if n_elements(event) eq 0 then begin
        rejected_index.add, candidate_id
    endif else begin
        events.add, event
        selected_index.add, candidate_id
    endelse
endforeach
selected_index = selected_index.toarray()
rejected_index = rejected_index.toarray()

root_dir = join_path([project.plot_dir,'diagnostic_plot'])
in_dir = join_path([root_dir,'all_subgroups'])
plot_files = file_search(join_path([in_dir,'*.pdf']))
plot_base_files = fgetbase(plot_files)
if n_elements(plot_files) ne candidates.length then stop

foreach the_dir, ['selected','rejected'] do begin
    the_index = (the_dir eq 'selected')? selected_index: rejected_index
    dirname = join_path([root_dir,'filter_df_group',the_dir])
    if file_test(dirname,/directory) eq 1 then file_delete, dirname, /recursive
    file_mkdir, dirname
    foreach candidate, candidates[the_index], id do begin
        file_suffix = 'azim_df_event_'+strjoin(time_string(candidate.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
        index = where(plot_base_files eq file_suffix, count)
        if count eq 0 then stop
        out_file = join_path([dirname,file_suffix])
        in_file = join_path([in_dir,file_suffix])
        file_copy, in_file, out_file
    endforeach
endforeach



;dirname = 'all_subgroups'
;azim_df_subgroup_gen_diagnostic_plot, events, dirname=dirname, project=project
stop

end
