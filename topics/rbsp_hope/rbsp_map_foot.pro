
pro rbsp_map_foot, tr, no_load_spice = no_load_spice, south = south

datetime = tr[0]
dur = (tr[1]-tr[0])/3600D
datetime0 = tr[0]-3600D
dur0 = dur+1.0

timespan, datetime0, dur0, /hours

get_timespan, t
print, time_string(t[0])
print, time_string(t[1])

omni_hro_load

get_data, 'OMNI_HRO_1min_BY_GSM', data = by
store_data, 'BY_GSM', data = by
get_data, 'OMNI_HRO_1min_BZ_GSM', data = bz
store_data, 'BZ_GSM', data = bz
get_data, 'OMNI_HRO_1min_flow_speed', data = flow
store_data, 'Vp_tvar', data = flow
get_data, 'OMNI_HRO_1min_proton_density', data = ppd
store_data, 'Np_tvar', data = ppd
get_data, 'OMNI_HRO_1min_SYM_H', data = symh
store_data, 'dst_tvar', data = symh

tdegap, 'BY_GSM', /overwrite
tdeflag, 'BY_GSM', 'linear', /overwrite
tdegap, 'BZ_GSM', /overwrite
tdeflag, 'BZ_GSM', 'linear', /overwrite
tdegap, 'Vp_tvar', /overwrite
tdeflag, 'Vp_tvar', 'linear', /overwrite
tdegap, 'Np_tvar', /overwrite
tdeflag, 'Np_tvar', 'linear', /overwrite
tdegap, 'dst_tvar', /overwrite
tdeflag, 'dst_tvar', 'linear', /overwrite

store_data, 'omni_imf', data = ['BY_GSM','BZ_GSM']
get_tsy_params, 'dst_tvar', 'omni_imf', 'Np_tvar', 'Vp_tvar', 't01', /speed, /imf_yz

; position.
if not keyword_set(no_load_spice) then begin
    rbsp_load_spice_kernels
endif

timespan, datetime, dur, /hours
get_timespan, t

rbsp_load_spice_state, probe = 'a', coord = 'gse', /no_spice_load
time_clip, 'rbspa_state_pos_gse', time_double(t[0]), time_double(t[1]), newname = 'pos_gse'

; set outname of foot point.
outname = 'ifoot_t01'

if not keyword_set(south) then begin
    ttrace2iono, 'pos_gse', external_model = 't01', par = 't01_par', in_coord = 'gse', out_coord = 'gse', $
        newname = outname, trace_var_name = 'trace_n', /km
endif else begin
    ttrace2iono, 'pos_gse', external_model = 't01', par = 't01_par', in_coord = 'gse', out_coord = 'gse', $
        newname = outname, trace_var_name = 'trace_s', /km, /south
endelse

cotrans, outname, outname+'_gei', /gse2gei
cotrans, outname+'_gei', outname+'_geo', /gei2geo

get_data, outname+'_geo', data = data
glat = atan(data.y[*,2]/sqrt(data.y[*,0]^2+data.y[*,1]^2))*180/!pi
glon = atan(data.y[*,1], data.y[*,0])*180/!pi

store_data, outname+'_glat_glon', data = {x:data.x, y:[[glat],[glon]]}, dlim = {colors:[2,4], labels:['GLAT','GLON']}

loadct, 39

if not keyword_set(south) then begin
    map_set, 90, 0, /stereographic, /continents, /grid, /isotropic, /noerase, $
        title = 'rbspa foot point from '+time_string(t[0])+' to '+time_string(t[1]), limit = [50,0,90,360]
    map_continents, /coasts
endif else begin
    map_set,-90, 0, /stereographic, /continents, /grid, /isotropic, /noerase, $
        title = 'rbspa foot point from '+time_string(t[0])+' to '+time_string(t[1]), limit = [15,0,-90,360]
    map_continents, /coasts
endelse

plots, glon, glat, color = 160

end

tr = time_double(['2013-04-14/07:00','2013-04-14/10:00'])
rbsp_map_foot, tr
end