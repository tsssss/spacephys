;+
; Bin data.
;-


pro pflux_survey_bin_data, var_type, probe=probe, bin_info_name=bin_info_name, $
    test_bin=test_bin

    if n_elements(bin_info_name) eq 0 then bin_info_name = 'spherical_spatial'
    bin_info = pflux_survey_load_bin_info(bin_info_name, probe=probe)

    bin_vars = bin_info.vars
    bin_dims = bin_info.bin_dims
    bin_mesh = bin_info.bin_mesh
    bin_mesh_range = bin_info.bin_mesh_range


    project = pflux_survey_load_project()
    data_file = join_path([project.data_dir,'binned_data',project.name+'_'+bin_info.name+'.cdf'])
    if file_test(data_file) eq 0 then begin
        foreach bin_var, bin_vars, ii do begin
            the_var = 'vertice_dim'+string(ii,format='(I0)')
            data = bin_info[bin_var].vertices
            cdf_save_var, the_var, value=data, filename=data_file
            settings = dictionary($
                'unit', bin_info[bin_var].unit, $
                'fieldnam', 'vertices on '+bin_var, $
                'display_type', 'scalar')
            cdf_save_setting, settings, varname=the_var, filename=data_file
        endforeach

        bin_vertices = list()
        foreach bin_var, bin_vars do bin_vertices.add, bin_info[bin_var].vertices
        mesh_vertices = mesh_grid(bin_vertices)
        foreach bin_var, bin_vars, ii do begin
            the_var = 'mesh_dim'+string(ii,format='(I0)')
            data = mesh_vertices[ii]
            cdf_save_var, the_var, value=data, filename=data_file, save_as_is=1
            settings = dictionary($
                'unit', bin_info[bin_var].unit, $
                'fieldnam', 'meshed vertices on '+bin_var, $
                'display_type', 'scalar')
            cdf_save_setting, settings, varname=the_var, filename=data_file
        endforeach
    endif


    ; Load data.
    pflux_survey_load_data, var_type, probe=probe
    prefix = 'rbsp'+probe+'_'
    the_var = prefix+var_type
    data = get_var_data(the_var, times=times)
    ndim = n_elements(data[0,*])
    




if n_elements(test_bin) ne 0 then test_single_bin = 1

if keyword_set(test_single_bin) then begin
    foreach bin_range, bin_mesh_range, bin_mesh_id do begin
        flag = 0
        foreach range, test_bin, range_id do begin
            flag += total(abs(range-bin_range[range_id]))
        endforeach
        if flag eq 0 then begin
            the_bin_mesh_id = bin_mesh_id
            break
        endif
    endforeach
    if n_elements(the_bin_mesh_id) eq 0 then stop
endif

;stop

    ; Bin data.
    data_list = list()
    foreach bin_index, bin_mesh, bin_mesh_id do begin
if keyword_set(test_single_bin) then if bin_mesh_id ne the_bin_mesh_id then begin
    data_list.add, !null
    continue
endif
        if n_elements(bin_index) eq 0 then begin
            data_list.add, !null
            continue
        endif else begin
            data_list.add, data[bin_index,*]
        endelse
    endforeach

    ; Get statistical results.
    data_types = ['mean','stddev']
    ndata_type = n_elements(data_types)

    data_info = dictionary()
    foreach data_type, data_types do begin
        all_data = list()
        for dim_id=0,ndim-1 do begin
            the_data = list()
            foreach data, data_list, data_id do begin
if keyword_set(test_single_bin) then if data_id ne the_bin_mesh_id then begin
    the_data.add, !values.f_nan
    continue
endif
                if n_elements(data) eq 0 then begin
                    the_val = !values.f_nan
                endif else begin
                    good_data = data[*,dim_id]
                    the_val = call_function(data_type, good_data, /nan)
                endelse
                the_data.add, the_val
            endforeach
            the_data = the_data.toarray()
            all_data.add, reform(the_data, bin_dims)
        endfor
        
        
        the_data = all_data.toarray()
        if ndim eq 1 then data_info[data_type] = reform(the_data) else data_info[data_type] = transpose(the_data,shift(indgen(ndim+1),-1))
    endforeach


    foreach data_type, data_types do begin
        cdf_save_var, the_var+'_'+data_type, value=data_info[data_type], filename=data_file, save_as_is=1
    endforeach


end

var_types = ['b0_sm', 'pf_sm_norm']
var_types = ['b0_sm','pf_sm_norm','ion_t','de_sm','dedot0_sm']
var_types = ['de_sm']
var_types = ['mlt','dis','mlat']
var_types = ['b0_sm']
probe = 'a'



; empty bin.
mlt_range = [-10]+[-1,1]*0.5
mlat_range = -20+[0,5]
dis_range = 4.0+[0,0.5]

; weird bin.
mlt_range = [8]+[-1,1]*0.5
mlat_range = -20+[0,5]
dis_range = 5.5+[0,0.5]
test_bin = list()
test_bin.add, mlt_range
test_bin.add, mlat_range
test_bin.add, dis_range

foreach var_type, var_types do pflux_survey_bin_data, var_type, probe=probe;, test_bin=test_bin
end
