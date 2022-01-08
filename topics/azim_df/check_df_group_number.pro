;+
; Check the number of DF groups.
;-

project = azim_df_load_project()
search_step = 'df_group'
candidates = list()
candidates.add, azim_df_search_post_midn_events(search_step=search_step), /extract
candidates.add, azim_df_search_pre_midn_events(search_step=search_step), /extract


times = list()
foreach candidate, candidates do times.add, candidate.time_range[0]
index = sort(times.toarray())
candidates = candidates[index]
foreach candidate, candidates, id do candidate.id = id


overlap_list = list()
current_list = list()
isolated_list = list()
event_list = list()
foreach the_group, candidates, id do begin
    str_id = string(id,format='(I0)')
    lprmsg, 'current candidate: '+str_id
    if current_list.length ne 0 then begin
        pre_group = current_list[-1]
        pre_time = pre_group.time_range
        the_time = the_group.time_range
        if min(the_time) ge max(pre_time) then begin
            lprmsg, 'no overlap with previous candidate, pop current_list'
            if current_list.length gt 1 then begin
                lprmsg, 'add current_list to overlap_list'
                overlap_list.add, current_list
                event_list.add, current_list
            endif else begin
                lprmsg, 'add current_list to isolated_list'
                isolated_list.add, current_list
                event_list.add, current_list
            endelse
            lprmsg, 'reset current_list'
            current_list = list()
        endif
    endif
    lprmsg, 'add '+str_id+' to current_list'
    current_list.add, the_group
    if the_group.id eq candidates.length-1 then begin
        if current_list.length gt 1 then begin
            lprmsg, 'add current_list to overlap_list'
            overlap_list.add, current_list
            event_list.add, current_list
        endif else begin
            lprmsg, 'add current_list to isolated_list'
            isolated_list.add, current_list
            event_list.add, current_list
        endelse
    endif
endforeach


; Resolve overlap list.
in_dir = join_path([project.plot_dir,'diagnostic_plot','df_group'])
foreach events, overlap_list, id do begin
    out_dir = join_path([project.plot_dir,'diagnostic_plot','overlap_df_group',string(id+1,format='(I0)')])
    if file_test(out_dir) eq 0 then file_mkdir, out_dir
    foreach the_event, events do begin
        base_name = 'azim_df_event_'+strjoin(time_string(the_event.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
        in_file = join_path([in_dir,base_name])
        out_file = join_path([out_dir,base_name])
        file_copy, in_file, out_file, /overwrite
    endforeach
endforeach
stop






search_step = 'uniq_df'
candidates = list()
candidates.add, azim_df_search_post_midn_events(search_step=search_step), /extract
candidates.add, azim_df_search_pre_midn_events(search_step=search_step), /extract

; Check UT-MLT.
min_r_square = 0.9

azim_list = list()
the_key = 'obs_mlt'
foreach event, candidates do begin
        xxs = list()
        yys = list()
        foreach df, event.df_list do begin
            xxs.add, df.obs_time
            yys.add, df[the_key]
        endforeach
        xxs = xxs.toarray()
        yys = yys.toarray()
        fit_result = linfit(xxs, yys, yfit=yfit)
        ss_res = total((yys-yfit)^2)
        ss_tot = total((yys-mean(yys))^2)
        r_square = 1-ss_res/ss_tot
        if r_square ge min_r_square then azim_list.add, event
endforeach

; Check UT-Rxy.
rxy_list = list()
the_key = 'obs_rxy'
foreach event, candidates do begin
        xxs = list()
        yys = list()
        foreach df, event.df_list do begin
            xxs.add, df.obs_time
            yys.add, df[the_key]
        endforeach
        xxs = xxs.toarray()
        yys = yys.toarray()
        fit_result = linfit(xxs, yys, yfit=yfit)
        ss_res = total((yys-yfit)^2)
        ss_tot = total((yys-mean(yys))^2)
        r_square = 1-ss_res/ss_tot
        if r_square ge min_r_square then rxy_list.add, event
endforeach


print, azim_list.length
print, rxy_list.length


stop
;---Copy plots.
in_dir = join_path([project.plot_dir,'diagnostic_plot','azim_df_uniq_df'])
out_dir = join_path([project.plot_dir,'diagnostic_plot','azim_df_uniq_df','azim_list'])
if file_test(out_dir) eq 0 then file_mkdir, out_dir
foreach the_event, azim_list do begin
    base_name = 'azim_df_event_'+strjoin(time_string(the_event.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
    in_file = join_path([in_dir,base_name])
    out_file = join_path([out_dir,base_name])
    file_copy, in_file, out_file
endforeach

in_dir = join_path([project.plot_dir,'diagnostic_plot','azim_df_uniq_df'])
out_dir = join_path([project.plot_dir,'diagnostic_plot','azim_df_uniq_df','rxy_list'])
if file_test(out_dir) eq 0 then file_mkdir, out_dir
foreach the_event, rxy_list do begin
    base_name = 'azim_df_event_'+strjoin(time_string(the_event.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
    in_file = join_path([in_dir,base_name])
    out_file = join_path([out_dir,base_name])
    file_copy, in_file, out_file
endforeach


end
