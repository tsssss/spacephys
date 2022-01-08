;+
; Generate daily-based overview plot.
;-

pro azim_df_gen_survey_plot_for_themis_good_triangle, project=project

    if n_elements(project) eq 0 then project = azim_df_load_project()
    log_file = join_path([project.data_dir,'azim_df_search_themis_good_triangle_list.txt'])
    lines = read_all_lines(log_file)
    ntime_range = n_elements(lines)
    time_ranges = dblarr(ntime_range,2)
    for ii=0, ntime_range-1 do begin
        tline = lines[ii]
        pos = strpos(tline, ' ')
        time_ranges[ii,0] = time_double(strmid(tline,0,pos))
        time_ranges[ii,1] = time_double(strmid(tline,pos+4))
    endfor

    secofday = 86400d
    date_ranges = time_ranges-(time_ranges mod secofday)
    date_ranges[*,1] += secofday

    index = uniq(date_ranges[*,0])
    nplot = n_elements(index)
    start_dates = reform(date_ranges[index,0])
    for ii=0, nplot-1 do begin
        start_date = start_dates[ii]
        index = where(date_ranges[*,0] eq start_date)
        end_date = max(date_ranges[index,1])
        time_range = [start_date,end_date]
        azim_df_gen_survey_plot, time_range, highlight_time_ranges=time_ranges[index,*]
    endfor


end

width = 80.*60/10
foreach prefix, 'th'+['a','b','c','d','e']+'_' do begin
    get_data, prefix+'db_tilt', times, data, limits=lim
    store_data, prefix+'theta', times, data-smooth(data, width, /nan, /edge_truncate), limits=lim
endforeach

end
