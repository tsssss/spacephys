;+
; Bin given var in the CDF file and given mapping_index.
;-

pro global_efield_bin_data, var, probe=mission_probe, bin_method=bin_method, project=project

;---Check input.
    if n_elements(var) eq 0 then message, 'No input var ...'
    if n_elements(project) eq 0 then project = global_efield_load_project()

;---Load mapping_info and data.
    ; Return mapping_info.
    mapping_index = global_efield_load_mapping(probe=mission_probe, bin_method=bin_method, project=project, $
        dim_vertices=dim_vertices, mesh_vertices=mesh_vertices)
    cell_dims = size(mapping_index,/dimensions)
    ncell = product(cell_dims)


;---Write to out_cdf.
    out_cdf = join_path([project.data_dir,'binned_data',mission_probe+'_'+bin_method+'_binned_data_v01.cdf'])
    if file_test(out_cdf) eq 0 then begin
        path = fgetpath(out_cdf)
        if file_test(path) eq 0 then file_mkdir, path
        cdfid = cdf_create(out_cdf)

        cell_ndim = n_elements(dim_vertices)
        bin_info = project[bin_method]
        dim_units = bin_info.units
        dim_labels = bin_info.labels

        count = lon64arr(ncell)
        for ii=0, ncell-1 do count[ii] = n_elements(*mapping_index[ii])
        count = reform(count, [cell_dims])

    ;---Save global settings.
        settings = dictionary($
            'mission_probe', mission_probe, $
            'coord_name', coord_name, $
            'coord_type', coord_type)
        cdf_save_setting, settings, filename=cdfid

    ;---Save dim_vertices.
        foreach data, dim_vertices, ii do begin
            tvar = 'vertice_dim'+sgnum2str(ii)
            cdf_save_var, tvar, value=data, filename=cdfid
            settings = dictionary($
                'fieldnam', 'vertices on '+dim_labels[ii], $
                'unit', dim_units[ii], $
                'display_type', 'scalar')
            cdf_save_setting, settings, filename=cdfid, varname=tvar
        endforeach


    ;---Save mesh_vertices.
        foreach data, mesh_vertices, ii do begin
            tvar = 'mesh_dim'+sgnum2str(ii)
            cdf_save_var, tvar, value=data, filename=cdfid, save_as_is=1
            settings = dictionary($
                'filednam', 'meshed vertices on '+dim_labels[ii], $
                'unit', dim_units[ii], $
                'display_type', 'scalar')
            cdf_save_setting, settings, filename=cdfid, varname=tvar
        endforeach

    ;---Save count.
        tvar = 'time_count'
        cdf_save_var, tvar, value=count, filename=cdfid, save_as_is=1
        settings = dictionary($
            'unit', 'min', $
            'display_type', 'scalar')
        cdf_save_setting, settings, filename=cdfid, varname=tvar
    endif else cdfid = cdf_open(out_cdf)


;---Load and bin data.
    ; Load var into tplot as <prefix>_var.
    global_efield_load_data, var, probe=mission_probe, project=project
    varname = project[mission_probe].prefix+var
    get_data, varname, times, data, limits=settings
    data_ndim = n_elements(data[0,*])
    cell_data = make_array([ncell,data_ndim], value=!values.f_nan)
    for ii=0, ncell-1 do begin
        if n_elements(*mapping_index[ii]) eq 0 then continue
        for jj=0, data_ndim-1 do begin
            the_data = data[*mapping_index[ii],jj]
            cell_data[ii,jj] = median(the_data)
        endfor
    endfor
    cell_data = reform(cell_data, [cell_dims,data_ndim])
    cdf_save_var, var, value=cell_data, filename=cdfid, save_as_is=1
    cdf_save_setting, dictionary(settings), filename=cdfid, varname=var

    cdf_close, cdfid


end

bin_method = 'r_sm_sph'
mission_probes = ['polar','rbsp'+letters('b'),'th'+letters('e')]
vars = ['ele_n','b_sm']

bin_method = 'r_sm_sph'
mission_probes = ['th'+letters('e')]
vars = ['u_sm']

bin_method = 'r_sm_sph'
mission_probes = ['c'+['1','2','3','4']]
vars = ['ele_n','b_sm']

bin_method = 'r_sm_sph'
mission_probes = ['polar','th'+letters('e'),'rbsp'+letters('b')]
vars = ['ele_t','ion_t']
vars = ['pmag']

bin_method = 'r_sm_sph'
mission_probes = 'c'+['1','2','3','4']
vars = ['b_sm','ele_n']

bin_method = 'r_sm_sph'
mission_probes = ['th'+letters('e')]
vars = ['e_sm']

bin_method = 'r_sm_sph'
mission_probes = 'rbsp'+['a','b']
vars = ['pf_sm_norm']

bin_method = 'r_sm_sph'
mission_probes = 'rbsp'+['a','b']
vars = ['ele_t','ion_t']

foreach mission_probe, mission_probes do begin
    foreach var, vars do begin
        global_efield_bin_data, var, probe=mission_probe, bin_method=bin_method
    endforeach
endforeach
end
