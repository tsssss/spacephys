;+
; Write given DF group to file.
;-
pro azim_df_subgroup_write_file, df_group, filename=file

    max_probe_length = 5
    tab = constant('4space')
    tformat = 'YYYY-MM-DD/hh:mm:ss'

    event_id = df_group.id
    msg = 'event_id'+string(event_id,format='(I8)')
    lprmsg, '', file
    lprmsg, msg, file
    msg = tab+$
        strjoin(time_string(df_group.time_range,tformat=tformat),' ')+tab+$
        df_group.region+tab+$
        df_group.search_name+tab+$
        strjoin(df_group.probes,',')
    lprmsg, msg, file

    df_list = df_group.df_list
    msg = 'df_list'+tab+string(df_list.length,format='(I0)')
    lprmsg, msg, file
    foreach df, df_list do begin
        azim_df_vertex_write, df, filename=file, prefix=tab
    endforeach

    edge_list = df_group.edge_list
    msg = 'edge_list'+tab+string(edge_list.length,format='(I0)')
    lprmsg, msg, file
    foreach edge, edge_list do begin
        azim_df_edge_write, edge, filename=file, prefix=tab
    endforeach

    triad_list = df_group.triad_list
    msg = 'triad_list'+tab+string(triad_list.length,format='(I0)')
    lprmsg, msg, file
    foreach triad, triad_list do begin
        azim_df_triad_write, triad, filename=file, prefix=tab
    endforeach
end
