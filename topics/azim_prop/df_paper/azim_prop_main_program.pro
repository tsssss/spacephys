
;---Load project.
    ;project = azim_prop_load_project()

;---Find storm times.
    ;if ~project.haskey('storm_times') then azim_prop_find_storm, project

;---Filter with spacecraft apogee MLT.
    ;if ~project.haskey('candidate_list') then azim_prop_filter_by_apogee, project


;---After the dataset is established.
    ; Init project
    ;project = azim_prop_load_project()

    ; Load data
    ;azim_prop_load_data, project, /reload

    ; Cross correlation to calculate time lag, ref time.
    azim_prop_calc_time_lag, project

    ; Calculate 3-spacecraft timing.
    azim_prop_calc_2d_vel, project

    ; Generate plots for cross correlation.
    azim_prop_plot_mlt_timing_only, project

    ; Generate plots for 3-spacecraft timing.
    azim_prop_plot_2d_timing, project

    ; Generate plots showing the timg lag.
    azim_prop_plot_time_lag, project

    ; Generate a plot for mlt vs angle.
    azim_prop_plot_angle, project

    ; Print the table in the paper.
    azim_prop_print_table, project

end
