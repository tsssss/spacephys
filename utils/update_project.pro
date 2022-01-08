;+
; Update the project dictionary in tplot and on disk.
;-

pro update_project, project
    store_data, project.var, 0, project
    tplot_save, project.var, filename=project.file
end