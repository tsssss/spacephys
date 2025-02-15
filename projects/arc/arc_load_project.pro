

function arc_load_project

;---Generic project info.
    project_id = 'arc'
    project_info = project_load(project_id)

;---Add events.
    tr = ['2016-09-02','2016-09-02/07:00']
    event = project_add_event(project_info, time_range=tr, update=1)

    project_info = project_update(project_info)
    return, project_info

end