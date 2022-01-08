;+
; A shortcut for updating the project.
;-

pro azim_prop_update_project, project

    store_data, project.var, 0, project
    tplot_save, project.var, filename=project.file

end
