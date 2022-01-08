;+
; Load mapping_hash, an internal var to save the mapping index from data to bins.
;-


;+
; Divide a bin in the mapping_hash into 2 bins for a given dimension.
;-
pro global_efield_divide_bin, mapping_hash, boundary=boundary, keys=keys, dim=dim, data=data

    if n_elements(mapping_hash) eq 0 then message, 'No mapping_hash ...'
    if ~isa(mapping_hash,'hash') then message, 'mapping_hash is not a hash ...'
    if n_elements(keys) eq 0 then keys = mapping_hash.keys()
    nkey = n_elements(keys)
    if nkey eq 0 then message, 'No key in mapping_hash ...'
    if n_elements(data) eq 0 then message, 'No data ...'
    ndim = (size(data,/dimensions))[1]
    nnew_bin = n_elements(boundary)/2

    foreach key, keys do begin
        dict = mapping_hash[key]
        current_boundary = dict.boundary
        current_data_index = dict.data_index
        current_data = reform(data[current_data_index,dim])
        mapping_hash.remove, key
        dim_keys = strarr(ndim)
        for ii=0, ndim-1 do dim_keys[ii] = '['+sgnum2str(current_boundary[0,ii])+','+sgnum2str(current_boundary[1,ii])+')'

        for ii=0, nnew_bin-1 do begin
            new_boundary = current_boundary
            new_boundary[*,dim] = boundary[*,ii]
            new_keys = dim_keys
            new_keys[dim] = '['+sgnum2str(new_boundary[0,dim])+','+sgnum2str(new_boundary[1,dim])+')'

            index = where(current_data ge boundary[0,ii] and current_data lt boundary[1,ii], new_count)
            if new_count eq 0 then begin
                lprmsg, 'No data in dim: '+sgnum2str(dim)+', range: '+new_keys[dim]+' ...'
                continue
            endif

            new_data_index = current_data_index[index]
            new_key = strjoin(new_keys, ',')
            mapping_hash[new_key] = dictionary($
                'boundary', new_boundary, $
                'data_index', new_data_index)
                ; Do not save count anymore, count is really a result of binning the time.
                ; Time should be treated equally as other physical quantities.
        endfor
    endforeach

end


pro global_efield_convert_mapping_hash_to_array, mapping_hash, filename=mapping_file, $
        dim_vertices=dim_vertices, mesh_vertices=mesh_vertices, mapping_index=mapping_index

;---Get the vertices on each dimension.
    keys = mapping_hash.keys()
    vector_ndim = n_elements((mapping_hash[keys[0]])['boundary'])/2

    ; The vertices on each dimension.
    dim_vertices = list()
    for ii=0, vector_ndim-1 do dim_vertices.add, list()
    foreach key, keys do begin
        tmp = mapping_hash[key]
        boundary = tmp['boundary']
        for ii=0, vector_ndim-1 do dim_vertices[ii].add, boundary[*,ii]
    endforeach
    for ii=0, vector_ndim-1 do begin
        current_vertices = dim_vertices[ii].toarray()
        current_vertices = sort_uniq(current_vertices)
        dim_vertices[ii] = current_vertices
    endfor

;---Get the 3D vertices, then get the mapping index from time series to cell.
    mesh_vertices = mesh_grid(dim_vertices)         ; The vertices of the meshed grid.
    mesh_dims = size(mesh_vertices[0],/dimensions)  ; cell vertices.
    cell_dims = mesh_dims-1                         ; cell center.
    cell_index_converter = lon64arr(vector_ndim)+1
    for ii=0, vector_ndim-2 do cell_index_converter[ii+1] = cell_index_converter[ii]*cell_dims[ii]
    mapping_index = ptrarr(cell_dims,/allocate_heap)
    foreach key, keys do begin
        tmp = mapping_hash[key]
        boundary = tmp['boundary']
        pos_of_cell = lonarr(vector_ndim)
        for ii=0, vector_ndim-1 do pos_of_cell[ii] = where(dim_vertices[ii] eq boundary[0,ii])
        pos_1d = total(pos_of_cell*cell_index_converter)
        *mapping_index[pos_1d] = tmp.data_index
    endforeach

    path = fgetpath(mapping_file)
    if file_test(path,/directory) eq 0 then file_mkdir, path
    if file_test(path,/directory) eq 0 then message, 'Invalid path: '+path+' ...'
    save, mapping_index, dim_vertices, mesh_vertices, filename=mapping_file

end



function global_efield_load_mapping, probe=mission_probe, bin_method=bin_method, project=project, $
    ; output
    dim_vertices=dim_vertices, mesh_vertices=mesh_vertices

;---Check input.
    if n_elements(bin_method) eq 0 then message, 'No bin_method ...'
    if n_elements(project) eq 0 then project = global_efield_load_project()

    mapping_file = join_path([project.data_dir,'mapping_data',project.name+'_'+mission_probe+'_'+bin_method+'.sav'])
    mapping_hash_file = join_path([project.data_dir,'mapping_data',project.name+'_'+mission_probe+'_'+bin_method+'_hash.sav'])

;---If mapping_index exists, read it and return.
    if file_test(mapping_file) eq 1 then begin
        ; Using restore directly seems to be unstable, try the obj and see how it works.
        ; Should have mapping_index, dim_vertices, mesh_vertices.
        ; dim_vertices is a list of ndim, the dimensions are [a,b,c].
        ; mapping_index is an array of pointer, in [a,b,c].
        savobj = obj_new('IDL_Savefile', mapping_file)
        vars_in_file = strlowcase(savobj->names())
        foreach var, ['mapping_index','dim_vertices','mesh_vertices'] do begin
            index = where(vars_in_file eq var, count)
            if count eq 0 then begin
                errmsg = handle_error(objid=savobj, 'File does not have var: '+var+' ...')
                return, !null
            endif
            savobj->restore, var
        endforeach
        obj_destroy, savobj
        if n_elements(mapping_index) gt 0 then return, mapping_index
    endif

;---If mapping_hash exits, use it to get mapping_index
; Save results and return.
    if file_test(mapping_hash_file) eq 1 then begin
        ; Using restore directly seems to be unstable, try the obj and see how it works.
        savobj = obj_new('IDL_Savefile', mapping_hash_file)
        save_var = 'mapping_hash'
        index = where(strlowcase(savobj->names()) eq save_var, count)
        if count ne 0 then savobj->restore, save_var
        obj_destroy, savobj
        if isa(mapping_hash, 'hash') then begin
            global_efield_convert_mapping_hash_to_array, mapping_hash, filename=mapping_file, $
                dim_vertices=dim_vertices, mesh_vertices=mesh_vertices, mapping_index=mapping_index
            return, mapping_index
        endif
    endif


;---Now we need to calculate mapping_hash and mapping_index
; Save both and return mapping_index.

;---Get bin settings.
    pos = strpos(bin_method,'_',/reverse_search)
    coord_type = strmid(bin_method,pos+1)
    orb_suffix = strmid(bin_method,0,pos)

    ; Get the default bin_size if it is not available.
    ; bin_info needs: ndim, sizes, ranges, units, names.
    the_key = 'bin_methods'
    if ~project.haskey(the_key) then project[the_key] = list()
    bin_methods = project[the_key]
    if bin_methods.where(bin_method) eq !null then begin
        bin_methods.add, bin_method
    endif

    vector_ndim = 3
    index = where(['car','sph','cyl'] eq coord_type, count)
    if count eq 0 then message, 'Unkown coord type ...'

    ; Bin size.
    r_bin_size = 0.5
    lat_bin_size = 5.
    lon_bin_size = 5.
    case coord_type of
        'car': bin_sizes = replicate(r_bin_size,vector_ndim)
        'sph': bin_sizes = [r_bin_size,lat_bin_size,lon_bin_size]
        'cyl': bin_sizes = [r_bin_size,lon_bin_size,r_bin_size]
    endcase

    ; Bin range.
    settings = project[mission_probe]
    r_range = settings.orb_r_range
    lat_range = [-90,90]
    lon_range = [-180,180]
    rxy_range = settings.orb_rxy_range
    x_range = settings.orb_x_range
    y_range = settings.orb_y_range
    z_range = settings.orb_z_range
    case coord_type of
        'car': bin_ranges = [[x_range],[y_range],[z_range]]     ; [2,ndim].
        'sph': bin_ranges = [[r_range],[lat_range],[lon_range]]
        'cyl': bin_ranges = [[rxy_range],[lon_range],[z_range]]
    endcase

    ; Labels.
    case coord_type of
        'car': bin_labels = ['x','y','z']
        'sph': bin_labels = ['r','lat','lon']
        'cyl': bin_labels = ['r','lon','z']
    endcase

    ; Units.
    case coord_type of
        'car': bin_units = replicate('Re',vector_ndim)
        'sph': bin_units = ['Re','deg','deg']
        'cyl': bin_units = ['Re','deg','Re']
    endcase

    bin_info = dictionary($
        'sizes', bin_sizes, $
        'ranges', bin_ranges, $
        'labels', bin_labels, $
        'units', bin_units)
    project[bin_method] = bin_info
    update_project, project

;---Read orbit data.
    lprmsg, 'Read orbit data and time from the CDF file ...'
    global_efield_load_data, orb_suffix, probe=mission_probe, project=project
    prefix = project[mission_probe].prefix
    orb_var = prefix+orb_suffix
    get_data, orb_var, times, orb_data
    ; Convert to sph or cyl.
    deg = 180d/!dpi
    case coord_type of
        'sph': begin
            dis = snorm(orb_data)
            lat = asin(orb_data[*,2]/dis)*deg
            lon = atan(orb_data[*,1],orb_data[*,0])*deg
            orb_data = [[dis],[lat],[lon]]
        end
        'cyl': begin
            rxy = snorm(orb_data[*,0:1])
            lon = atan(orb_data[*,1],orb_data[*,0])*deg
            orb_data = [[rxy],[lon],[orb_data[*,2]]]
        end
        'car': ; do nothing.
    endcase

;---Make the bins on each dimension and the 3D mesh grid.
    lprmsg, 'Get the bins on each dimension and mesh them to 3D grids ...'
    dim_vertices = list()
    bin_ranges = bin_info.ranges
    bin_sizes = bin_info.sizes
    for ii=0, vector_ndim-1 do dim_vertices.add, make_bins(bin_ranges[*,ii],bin_sizes[ii])


;---Calculate and save mapping_hash.
    lprmsg, 'Calculate mapping from time series to bins, using: '+orb_var+' ...'
    current_boundary = fltarr(2,vector_ndim)
    for ii=0, vector_ndim-1 do current_boundary[*,ii] = minmax(dim_vertices[ii])
    current_keys = strarr(vector_ndim)
    for jj=0, vector_ndim-1 do current_keys[jj] = '['+sgnum2str(current_boundary[0,jj])+','+sgnum2str(current_boundary[1,jj])+')'
    current_key = strjoin(current_keys, ',')
    ntime = n_elements(times)
    mapping_hash = orderedhash()
    mapping_hash[current_key] = dictionary($
        'boundary', current_boundary, $
        'data_index', indgen(ntime,/l64))

    ; Bin data according to ranges in each dimension.
    for ii=0, vector_ndim-1 do begin
        ; In [2,nbin].
        current_boundary = dim_vertices[ii]
        current_boundary = transpose([[current_boundary[0:-2]],[current_boundary[1:-1]]])
        global_efield_divide_bin, mapping_hash, boundary=current_boundary, dim=ii, data=orb_data
    endfor

    ; Save mapping_hash.
    lprmsg, 'Save mapping_hash to '+mapping_hash_file+' ...'
    path = fgetpath(mapping_hash_file)
    if file_test(path,/directory) eq 0 then file_mkdir, path
    if file_test(path,/directory) eq 0 then message, 'Invalid path: '+path+' ...'
    save, mapping_hash, filename=mapping_hash_file


; Calculate and save mapping_index.
    lprmsg, 'Convert mapping_hash, and save it to '+mapping_file+' ...'
    global_efield_convert_mapping_hash_to_array, mapping_hash, filename=mapping_file, $
        dim_vertices=dim_vertices, mesh_vertices=mesh_vertices, mapping_index=mapping_index
    return, mapping_index


end

mission_probes = ['polar','rbsp'+letters('b'),'th'+letters('e')]
mission_probes = ['th'+letters(['d','e'])]
bin_methods = 'r_sm_'+['car','sph']
foreach mission_probe, mission_probes do begin
    foreach bin_method, bin_methods do begin
        mapping_index = global_efield_load_mapping(probe=mission_probe, bin_method=bin_method)
    endforeach
endforeach
end
