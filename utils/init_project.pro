;+
; Initialize a project dictionary using a project name.
; This is a standard dictionary includes
;   name. A string of the project name.
;   var. A string of tplot var for the project dictionary.
;   root_dir. A string of the project's root directory.
;   plot_dir. A string of path to save plots.
;   data_dir. A string of path to save data.
;   code_dir. A string of path to the programs.
;   file. A string of full file name to save the project dictionary.
;   constant. A dictionary of commonly used constants (deg,rad,re,rgb,xyz,secofday)
;
; Note: the program does not save the project dictionary to file.
;
; project_name. A string of project name.
;-

function init_project, project_name, errmsg=errmsg, no_save=no_save
    errmsg = ''
    retval = !null

    if n_elements(project_name) eq 0 then begin
        errmsg = handle_error('No project name ...')
        return, retval
    endif

    project = dictionary()
    project['name'] = strlowcase(project_name[0])
    project['var'] = project.name+'_project_info'
    project['root_dir'] = join_path([googledir(),'works',project.name])
    project['data_dir'] = join_path([project.root_dir,'data'])
    project['plot_dir'] = join_path([project.root_dir,'plot'])
    project['code_dir'] = join_path([homedir(),'Projects','idl','spacephys','topics',project.name])
    project['file'] = join_path([project.data_dir,project.name+'_project_info.tplot'])

    if keyword_set(no_save) then return, project


    lprmsg, 'Initializing project: '+project.name+' ...'
    lprmsg, 'Check root directory ...'
    if file_test(project.root_dir,/directory) eq 0 then begin
        lprmsg, 'Creating root directory ...'
        file_mkdir, project.root_dir
        if file_test(project.root_dir,/directory) eq 0 then begin
            msg = handle_error('Fail to create root directory ...')
            lprmsg, msg
            return, retval
        endif else lprmsg, 'Root directory created ...'
    endif else lprmsg, 'Root directory exists ...'
    if file_test(project.data_dir,/directory) eq 0 then begin
        lprmsg, 'Creating data directory ...'
        file_mkdir, project.data_dir
    endif else lprmsg, 'Data directory exists ...'
    if file_test(project.plot_dir,/directory) eq 0 then begin
        lprmsg, 'Creating plot directory ...'
        file_mkdir, project.plot_dir
    endif else lprmsg, 'Plot directory exists ...'
    if file_test(project.code_dir,/directory) eq 0 then begin
        lprmsg, 'Creating code directory ...'
        file_mkdir, project.code_dir
    endif else lprmsg, 'Code directory exists ...'

    update_project, project

    lprmsg, 'The project dictionary is saved to file ...'
    return, project

end

project = init_project('global_efield')
end
