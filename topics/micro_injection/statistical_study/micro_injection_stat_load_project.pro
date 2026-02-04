function micro_injection_stat_load_project

    project_id = 'micro_injection'
    project_info = project_load(project_id)

    ; sample_energys = [48,104d]      ; keV. Works for both MMS/FEEPS and THEMIS/SST.
    ; 48 is too clost to the lowest energy bin of MMS/FEEPS.
    sample_energys = [67,122d]      ; keV. Works for both MMS/FEEPS and THEMIS/SST.
    project_info['sample_energys'] = sample_energys
    project_info['local_root'] = join_path([default_local_root(),'sdata','micro_injection'])
    project_info['data_dir'] = join_path([default_local_root(),'sdata','micro_injection'])

    scale_info = dictionary({s0:40d, s1:4000, dj:1d/8, ns:0d })
    project_info['scale_info'] = scale_info
    project_info['freq_range'] = minmax(1d/([1d,10]*60))  ; in Hz.
    project_info['common_time_step'] = 20d  ; sec.
    project_info['pdyn_range'] = [2d,10]

    return, project_info

end