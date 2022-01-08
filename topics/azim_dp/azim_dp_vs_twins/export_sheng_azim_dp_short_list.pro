
;---Get all events, from fig_v2d.
    out_file = join_path([srootdir(),'sheng_azim_dp_short_list.txt'])
    if file_test(out_file) eq 1 then file_delete, out_file
    ftouch, out_file
    
    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_find_dfgroup(project=project)
    short_list = list()
    foreach event, events do begin
        direction = (strsplit(event.region,'%',/extract))[1]

        if direction ne 'eastward' and direction ne 'westward' then continue
        short_list.add, event
    endforeach

    foreach event, short_list, event_id do begin
        event.id = event_id+1
        azim_df_subgroup_write_file, event, filename=out_file
    endforeach

    stop
end
