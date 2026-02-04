
function micro_injection_load_project

;---Generic project info.
    project_id = 'micro_injection'
    project_info = project_load(project_id)

;---Add events.
    tr = ['2015-09-05/14:00','2015-09-05/18:00']
    event = project_add_event(project_info, time_range=tr)
    event['pad_times'] = time_double([$
        '2015-09-05/11:38:20' ])

    tr = ['2015-09-11/18:30','2015-09-11/20:40']
    event = project_add_event(project_info, time_range=tr)
    event['pad_times'] = time_double([$
        '2015-09-01/11:38:20' ])

    tr = ['2015-09-01/11:30','2015-09-01/13:10']
    event = project_add_event(project_info, time_range=tr)
    event['pad_times'] = time_double([$
        '2015-09-01/11:38:20', $
        '2015-09-01/11:41:50', $
        '2015-09-01/11:44:10', $
        '2015-09-01/11:49:10', $
        '2015-09-01/11:56:15', $
        '2015-09-01/11:58:25', $
        '2015-09-01/12:00:50', $
        '2015-09-01/12:03:50', $
        '2015-09-01/12:06:10', $
        '2015-09-01/12:08:30', $
        '2015-09-01/12:12:25', $
        '2015-09-01/12:15:00', $
        '2015-09-01/12:29:40', $
        '2015-09-01/12:34:50', $
        '2015-09-01/12:39:40', $
        '2015-09-01/12:46:05', $
        '2015-09-01/12:50:55', $
        '2015-09-01/12:56:10' ])
    
    tr = ['2015-09-01/18:00','2015-09-01/21:30']
    event = project_add_event(project_info, time_range=tr)
    event['pad_times'] = time_double([$
        '2015-09-01/18:12:40', $
        '2015-09-01/18:18:00', $
        '2015-09-01/19:15:20', $
        '2015-09-01/19:03:20', $
        '2015-09-01/18:56:40', $
        '2015-09-01/19:19:40', $
        '2015-09-01/19:24:20', $
        '2015-09-01/19:29:40', $
        '2015-09-01/19:36:40' ])
    
    tr = ['2015-09-01/10:00','2015-09-01/23:30']
    event = project_add_event(project_info, time_range=tr)
    event['pad_times'] = time_double([$
        '2015-09-01/18:12:40', $
        '2015-09-01/19:36:40' ])
    
    tr = ['2016-03-05/17:30','2016-03-05/21:30']
    event = project_add_event(project_info, time_range=tr)
    event['pad_times'] = time_double([$
        '2016-03-05/18:12:40', $
        '2016-03-05/18:13:40' ])


    tr = ['2017-01-12/00:00','2017-01-12/03:00']
    event = project_add_event(project_info, time_range=tr)
    event['pad_times'] = time_double([$
        '2017-01-12/01:12:40', $
        '2017-01-12/01:13:40' ])



    project_info = project_update(project_info)
    return, project_info

end