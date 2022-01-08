;+
; Return a dictionary of settings of the project.
;-

function global_efield_load_project

    project = load_project('global_efield')
    project_updated = 0

;---Check overall settings.
    the_key = 'time_step'
    if ~project.haskey(the_key) then begin
        project[the_key] = 60d ; 1 min.
        project_updated = 1
    endif

;---Check mission based settings.
    mission_probes = ['polar','rbsp'+letters('b'),'th'+letters('e'),'c'+['1','2','3','4']]
    the_key = 'all_mission_probes'
    if ~project.haskey(the_key) then begin
        project[the_key] = mission_probes
        project_updated = 1
    endif

    foreach mission_probe, mission_probes do begin
        if project.haskey(mission_probe) then continue
        project_updated = 1
        settings = dictionary()
        case mission_probe of
            'polar': begin
                settings['prefix'] = 'po_'
                settings['routine_name'] = 'polar'
                settings['probe'] = ''
                settings['year_range'] = [1996,2008]
                settings['time0'] = '1996-02-28'
                settings['orb_x_max'] = 11.
                settings['orb_y_max'] = 11.
                settings['orb_z_max'] = 11.
            end
            'rbspa': begin
                settings['prefix'] = 'rbspa_'
                settings['routine_name'] = 'rbsp'
                settings['probe'] = 'a'
                settings['year_range'] = [2013,2018]
                settings['time0'] = '2012-09-05'
                settings['orb_x_max'] = 7.
                settings['orb_y_max'] = 7.
                settings['orb_z_max'] = 5.
            end
            'rbspb': begin
                settings['prefix'] = 'rbspb_'
                settings['routine_name'] = 'rbsp'
                settings['probe'] = 'b'
                settings['year_range'] = [2013,2018]
                settings['time0'] = '2012-09-05'
                settings['orb_x_max'] = 7.
                settings['orb_y_max'] = 7.
                settings['orb_z_max'] = 5.
            end
            'tha': begin
                settings['prefix'] = 'tha_'
                settings['routine_name'] = 'themis'
                settings['probe'] = 'a'
                settings['year_range'] = [2008,2018]
                settings['time0'] = '2008-01-01'
                settings['orb_x_max'] = 16.
                settings['orb_y_max'] = 16.
                settings['orb_z_max'] = 10.
            end
            'thb': begin
                settings['prefix'] = 'thb_'
                settings['routine_name'] = 'themis'
                settings['probe'] = 'b'
                settings['year_range'] = [2008,2018]
                settings['time0'] = '2008-01-01'
                settings['orb_x_max'] = 40.
                settings['orb_y_max'] = 40.
                settings['orb_z_max'] = 20.
            end
            'thc': begin
                settings['prefix'] = 'thc_'
                settings['routine_name'] = 'themis'
                settings['probe'] = 'c'
                settings['year_range'] = [2008,2018]
                settings['time0'] = '2008-01-01'
                settings['orb_x_max'] = 30.
                settings['orb_y_max'] = 30.
                settings['orb_z_max'] = 15.
            end
            'thd': begin
                settings['prefix'] = 'thd_'
                settings['routine_name'] = 'themis'
                settings['probe'] = 'd'
                settings['year_range'] = [2008,2018]
                settings['time0'] = '2008-01-01'
                settings['orb_x_max'] = 16.
                settings['orb_y_max'] = 16.
                settings['orb_z_max'] = 10.
            end
            'the': begin
                settings['prefix'] = 'the_'
                settings['routine_name'] = 'themis'
                settings['probe'] = 'e'
                settings['year_range'] = [2008,2018]
                settings['time0'] = '2008-01-01'
                settings['orb_x_max'] = 16.
                settings['orb_y_max'] = 16.
                settings['orb_z_max'] = 10.
            end
            'c1': begin
                settings['prefix'] = 'c1_'
                settings['routine_name'] = 'cluster'
                settings['probe'] = '1'
                settings['year_range'] = [2001,2018]
                settings['time0'] = '2001-01-30'
                settings['orb_x_max'] = 20.
                settings['orb_y_max'] = 20.
                settings['orb_z_max'] = 20.
            end
            'c2': begin
                settings['prefix'] = 'c2_'
                settings['routine_name'] = 'cluster'
                settings['probe'] = '2'
                settings['year_range'] = [2001,2018]
                settings['time0'] = '2001-01-30'
                settings['orb_x_max'] = 20.
                settings['orb_y_max'] = 20.
                settings['orb_z_max'] = 20.
            end
            'c3': begin
                settings['prefix'] = 'c3_'
                settings['routine_name'] = 'cluster'
                settings['probe'] = '3'
                settings['year_range'] = [2001,2018]
                settings['time0'] = '2001-01-30'
                settings['orb_x_max'] = 20.
                settings['orb_y_max'] = 20.
                settings['orb_z_max'] = 20.
            end
            'c4': begin
                settings['prefix'] = 'c4_'
                settings['routine_name'] = 'cluster'
                settings['probe'] = '4'
                settings['year_range'] = [2001,2018]
                settings['time0'] = '2001-01-30'
                settings['orb_x_max'] = 20.
                settings['orb_y_max'] = 20.
                settings['orb_z_max'] = 20.
            end
        endcase
        settings['years'] = string(make_bins(settings['year_range'],1),format='(I4)')
        settings['time0'] = time_double(settings['time0'])
        settings['orb_rxy_max'] = max([settings['orb_x_max'],settings['orb_y_max']])
        settings['orb_r_max'] = max([settings['orb_rxy_max'],settings['orb_z_max']])
        settings['orb_x_range'] = [-1,1]*settings['orb_x_max']
        settings['orb_y_range'] = [-1,1]*settings['orb_y_max']
        settings['orb_z_range'] = [-1,1]*settings['orb_z_max']
        settings['orb_r_range'] = [1,settings['orb_r_max']]
        settings['orb_rxy_range'] = [1,settings['orb_rxy_max']]

        project[mission_probe] = settings
    endforeach

    if project_updated then update_project, project
    return, project

end
