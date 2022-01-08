;+
; Read DF group.
;-

function azim_df_subgroup_read_file, file
test = 0

    lines = read_all_lines(file)
    nline = n_elements(lines)

    df_group_list = list()
    for line_id=0,nline-1 do begin
        tline = lines[line_id]
        if keyword_set(test) then lprmsg, tline
        if strtrim(tline) eq '' then continue
        infos = strsplit(tline, ' ', /extract)
        ninfo = n_elements(infos)

        if infos[0] eq 'event_id' then begin
            df_group = dictionary('id',fix(infos[1]))
            line_id += 1
            tline = lines[line_id]
            infos = strsplit(tline, ' ', /extract)
            df_group['time_range'] = time_double(infos[0:1])
            df_group['region'] = infos[2]
            df_group['search_name'] = infos[3]
            df_group['probes'] = strsplit(infos[4],',',/extract)
            df_group['df_list'] = list()
            df_group['triad_list'] = list()
            df_group['edge_list'] = list()
            continue
        endif

        if infos[0] eq 'df_list' then begin
            ndf = fix(infos[1])
            df_list = df_group.df_list
            for ii=0,ndf-1 do begin
                line_id += 1
                df_list.add, azim_df_vertex_read(lines[line_id])
            endfor
            continue
        endif

        if infos[0] eq 'edge_list' then begin
            nedge = fix(infos[1])
            edge_list = df_group.edge_list
            for ii=0,nedge-1 do begin
                line_id += 1
                edge_list.add, azim_df_edge_read(lines[line_id])
            endfor
            continue
        endif

        if infos[0] eq 'triad_list' then begin
            ntriad = fix(infos[1])
            triad_list = df_group.triad_list
            for ii=0,ntriad-1 do begin
                line_id += 1
                triad_list.add, azim_df_triad_read(lines[line_id])
            endfor
        endif

        df_group_list.add, df_group
    endfor

    return, df_group_list
end


fn = '/Volumes/GoogleDrive/My Drive/works/works/azim_df/data/azim_df_search_event_pre_midn_beyond_15Re_search_df_group.txt'
events = azim_df_subgroup_read_file(fn)
end
