;+
; Return a project dictionary.
;
; project_info. A string for the project name, or the project dictionary.
;   The project related initialization function must exists, otherwise
;   the program will fail.
;-
function load_project, project_info
    errmsg = ''
    retval = !null


    case strlowcase(typename(project_info)) of
        'string': project = init_project(strlowcase(project_info[0]),/no_save)
        'dictionary': project = project_info
        else: errmsg = 'Invalid type of project_info ...'
    endcase

    if errmsg eq '' then if ~project.haskey('file') then $
        errmsg = 'Invalid project dictionary ...'
    if errmsg ne '' then begin
        errmsg = handle_error(errmsg)
        return, retval
    endif
    
    if file_test(project.file) eq 0 then update_project, project
    tplot_restore, filename=project.file
    project = get_var_data(project.var)
    
    ; GoogleDrive appears as different names in different OS.
    if ~file_test(project.root_dir,/directory) then begin
        project.root_dir = join_path([googledir(),'works','works',project.name])
        project.data_dir = join_path([project.root_dir,'data'])
        project.plot_dir = join_path([project.root_dir,'plot'])
        project.code_dir = join_path([homedir(),'Project','idl','spacephys','topics',project.name])
        project.file = join_path([project.data_dir,project.name+'_project_info.tplot'])
    endif
    if ~file_test(project.file) then call_procedure, project.name+'_init_project', project
    if ~file_test(project.file) then message, 'Something wrong in loading the project info ...'

    return, project
end
