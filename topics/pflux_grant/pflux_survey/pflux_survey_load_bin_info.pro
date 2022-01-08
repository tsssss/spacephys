;+
; Create mapping_index for given vars. Adopted from global_efield_load_mapping.
;
; bin_info should be a dict of {var, ranges}.
;-

function pflux_survey_load_bin_info_mesh_per_dim, old_list, bins

    new_list = list()
    foreach the_list, old_list do begin
        foreach bin, bins do begin
            ;lprmsg, '    '+strjoin(string(bin.range,format='(F5.1)'),' to ')+' ...'
            new_list.add, set_inter(the_list, bin.index)
        endforeach
    endforeach

    return, new_list

end


function pflux_survey_load_bin_info_mesh_per_dim_for_range, old_list, bins

    new_list = list()
    foreach the_list, old_list do begin
        foreach bin, bins do begin
            ;lprmsg, '    '+strjoin(string(bin.range,format='(F5.1)'),' to ')+' ...'
            tmp_list = list(bin.range)
            tmp_list.add, the_list, /extract
            new_list.add, tmp_list
        endforeach
    endforeach

    return, new_list

end


function pflux_survey_load_bin_info, bin_info, probe=probe, project=project, $
    ; output
    dim_vertices=dim_vertices, mesh_vertices=mesh_vertices

    if n_elements(project) eq 0 then project = pflux_survey_load_project()
    prefix = 'rbsp'+probe+'_'
    if size(bin_info,/type) eq 7 then begin
        bin_info_name = bin_info
    endif else begin
        bin_info_name = bin_info.name
    endelse

    bin_info_file = join_path([project.data_dir,'mapping_data',project.name+'_'+prefix+'bin_info_'+bin_info_name+'.sav'])
    path = fgetpath(bin_info_file)
    if file_test(path,/directory) eq 0 then file_mkdir, path

;---If mapping_index exists, read it and return.
    if file_test(bin_info_file) eq 1 then begin
        ; Using restore directly seems to be unstable, try the obj and see how it works.
        ; Should have mapping_index, dim_vertices, mesh_vertices.
        ; dim_vertices is a list of ndim, the dimensions are [a,b,c].
        ; mapping_index is an array of pointer, in [a,b,c].
        savobj = obj_new('IDL_Savefile', bin_info_file)
        vars_in_file = strlowcase(savobj->names())
        foreach var, ['bin_info'] do begin
            index = where(vars_in_file eq var, count)
            if count eq 0 then begin
                errmsg = handle_error(objid=savobj, 'File does not have var: '+var+' ...')
                return, !null
            endif
            savobj->restore, var
        endforeach
        obj_destroy, savobj
        return, bin_info
    endif


    if size(bin_info,/type) eq 7 then message, 'No input bin_info ...'

;---Get bin settings.
    ; Load common time.
    time_var = 'time'
    data_file = join_path([project.data_dir,project.data_file_suffix])
    times = cdf_read_var(time_var, filename=data_file)
    ntime = n_elements(times)

    vars = bin_info.vars
    ndim = n_elements(vars)
    bin_dims = list()
    dim_vertices = list()
    foreach var, vars do begin
        lprmsg, 'Binning '+var+' ...'
        ; Load data.
        pflux_survey_load_data, var, probe=probe
        the_var = prefix+var
        data = get_var_data(the_var)
        if n_elements(data) ne ntime then message, 'Invalid bin_var ...'

        ; Bin data.
        the_info = bin_info[var]
        nbin = n_elements(the_info.vertices)-1
        bin_dims.add, nbin
        bins = list()
        for ii=0, nbin-1 do begin
            range = the_info.vertices[ii:ii+1]
            relation = '[)'
            if range[1] lt range[0] then begin
                relation = ')['
                range = range[[1,0]]
            endif
            index = lazy_where(data, relation, range, count=count)
            if count eq 0 then index = !null
            bins.add, dictionary($
                'range', range, $
                'count', count, $
                'index', index)
        endfor
        the_info['bins'] = bins
    endforeach
    bin_dims = bin_dims.toarray()
    bin_info['bin_dims'] = bin_dims


;---Mesh bins.
    bin_mesh = list()
    bin_mesh.add, findgen(ntime)
    bin_mesh_range = list()
    bin_mesh_range.add, list()
    foreach var, reverse(vars) do begin ; Need to reverse to make reform work for 3d.
        lprmsg, 'Meshing '+var+' ...'
        bin_mesh = pflux_survey_load_bin_info_mesh_per_dim(bin_mesh, bin_info[var].bins)
        bin_mesh_range = pflux_survey_load_bin_info_mesh_per_dim_for_range(bin_mesh_range, bin_info[var].bins)
    endforeach
    bin_info['bin_mesh'] = bin_mesh
    bin_info['bin_mesh_range'] = bin_mesh_range

    ; Counts per meshed bin.
    counts = list()
    foreach data, bin_mesh do counts.add, n_elements(data)
    counts = counts.toarray()
    counts_3d = reform(counts, bin_dims)
    bin_info['bin_mesh_count'] = counts_3d


;---Save data.
    save, bin_info, filename=bin_info_file
    return, bin_info

end

mlt_vertices = [smkarthm(-11.5,11.5,1, 'dx'),-11.5]
bin_info = dictionary($
    'name', 'spherical_spatial', $
    'vars', ['mlt','mlat','dis'], $
    'mlt', dictionary($
        'vertices', mlt_vertices, $
        'unit', 'hr'), $
    'mlat', dictionary($
        'vertices', smkarthm(-20,20,5, 'dx'), $
        'unit', 'deg'), $
    'dis', dictionary($
        'vertices', smkarthm(4,6,0.5, 'dx'), $
        'unit', 'Re'))

probe = 'b'
tmp = pflux_survey_load_bin_info(bin_info, probe=probe)

end