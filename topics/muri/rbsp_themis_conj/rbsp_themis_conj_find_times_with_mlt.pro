;+
; Load Spacecraft positions and find conjunctions based on MLT.
;-

function find_apogee, r_var, get_name=get_name

    out_var = r_var+'_apogee'
    if keyword_set(get_name) then return, out_var

    get_data, r_var, times, r_coord, limits=lim
    dis = snorm(r_coord)
    max_dis = max(dis)
    index = where(dis ge max_dis*0.5, count)
    if count eq 0 then return, ''

    sectors = time_to_range(index,time_step=1)
    nsector = n_elements(sectors)*0.5
    apogee_index = dblarr(nsector)
    for ii=0,nsector-1 do begin
        the_index = make_bins(sectors[ii,*],1)
        tmp = max(dis[the_index], index)
        apogee_index[ii] = the_index[index]
    endfor
    apogee_index = sort_uniq(apogee_index)

    store_data, out_var, times[apogee_index], r_coord[apogee_index,*], limits=lim
    return, out_var
    
end

rbsp_probes = ['a','b']
themis_probes = ['a','d','e']
time_range = ['2013-01-01','2019-01-01']
coord = 'sm'
orbit_vars = list()
foreach probe, rbsp_probes do begin
    prefix = 'rbsp'+probe+'_'
    var = prefix+'r_'+coord
    apogee_var = find_apogee(var, get_name=1)
    if check_if_update(apogee_var) then begin
        var = ml_rbsp_read_pos(time_range, probe=probe, coord=coord)
        apogee_var = find_apogee(var)
    endif
    get_data, apogee_var, times, r_sm
    index = uniq(times, sort(times))
    times = times[index]
    r_sm = r_sm[index,*]
    mlt = pseudo_mlt(r_sm)
    store_data, prefix+'mlt', times, mlt
endforeach

foreach probe, themis_probes do begin
    prefix = 'th'+probe+'_'
    var = prefix+'r_'+coord
    apogee_var = find_apogee(var, get_name=1)
    if check_if_update(apogee_var) then begin
        var = themis_read_orbit(time_range, probe=probe, coord=coord)
        apogee_var = find_apogee(var)
    endif
    get_data, apogee_var, times, r_sm
    index = uniq(times, sort(times))
    times = times[index]
    r_sm = r_sm[index,*]
    mlt = pseudo_mlt(r_sm)
    store_data, prefix+'mlt', times, mlt
endforeach


probes = ['rbsp'+rbsp_probes,'th'+themis_probes]
colors = sgcolor(['blue','cyan','red','green','purple'])
vars = probes+'_mlt'
var = 'mlt_combo'
stplot_merge, vars, newname=var, labels=strupcase(probes), colors=colors
options, var, 'labflag', -1
options, var, 'yrange', [-1,1]*12
options, var, 'ystyle', 1
options, var, 'constant', [0,[-1,1]*3,[-1,1]*6]
options, var, 'ytitle', 'Apogee MLT (h)'

end