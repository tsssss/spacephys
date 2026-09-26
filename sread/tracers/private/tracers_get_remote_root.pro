
;+
; :Purpose: Get remote root, username, and password.
; :Returns: Remote info as requested.
; :Arguments:
;   id: in, optional, 'username' or 'password'
;-
function tracers_get_remote_root, id
    compile_opt idl2

    if n_elements(id) eq 0 then id = 'url'

    usrname = 'tracers-sot'
    pwd = 'SciOpsTeamFlight!'
    if id eq 'username' then return, usrname
    if id eq 'password' then return, pwd
    
    url = 'https://tracers-portal.physics.uiowa.edu/teams/flight/'
    return, url

end
