;+
; Load project info.
;-

function rbsp_themis_injection_load_project, reset=reset

    project_name = 'rbsp_themis_injection'
    project = load_project(project_name)
    project_updated = 0

    if keyword_set(reset) then begin
        del_data, project.var
        file_delete, project.file, /allow_nonexistent
        project = load_project(project_name)
    endif

;---Overall settings.
    the_key = 'search_time_range'
    if ~project.haskey(the_key) then begin
        search_time_range = time_double(['2013','2018'])
        project[the_key] = search_time_range
        project_updated = 1
    endif

    if project_updated then update_project, project
    return, project
end


project_info = rbsp_themis_injection_load_project()
end