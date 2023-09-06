;+
; Show storm dependence.
;-

test = 0


probes = ['a','b']
time_range = time_double(['2012-10','2015-10'])
foreach probe, probes do begin
    prefix = 'rbsp'+probe+'_'
    
    the_var = prefix+'orbit_times'
    if check_if_update(the_var) then begin
        orbit_time_ranges = pflux_grant_calculate_orbit_time_range(time_range, probe=probe)
        store_data, the_var, 0, orbit_time_ranges
    endif
    orbit_time_ranges = get_var_data(the_var)
    orbit_times = (orbit_time_ranges[*,0]+orbit_time_ranges[*,1])*0.5
    norbit = n_elements(orbit_times)



    pflux_survey_load_data, 'pf_fac_norm', probe=probe
    get_data, prefix+'pf_fac_norm', times, pf

    min_dis = 4.
    pflux_survey_load_data, 'dis', probe=probe
    get_data, prefix+'dis', times, dis
    index = where(dis le min_dis)
    pf[index,*] = !values.f_nan
    pf = snorm(pf)

    pflux_survey_load_data, 'mlt', probe=probe
    apogee_mlts = get_var_data(prefix+'mlt', at=orbit_times)
    the_var = prefix+'apogee_mlt'
    store_data, the_var, orbit_times, apogee_mlts
    add_setting, the_var, /smart, dictionary($
        'display_type', 'scalar', $
        'unit', 'hr', $
        'short_name', 'MLT')

    data = fltarr(norbit)
    for ii=0,norbit-1 do begin
        index = where_pro(times, '[)', orbit_time_ranges[ii,*], count=count)
        duration = total(orbit_time_ranges[ii,*]*[-1,1])
        if count eq 0 then begin
            data[ii] = !values.f_nan
        endif else begin
            data[ii] = mean(pf[index,*], /nan);*duration
        endelse
    endfor

    ; Adjust unit.
    unit = 'mW/m!U2!N'
    
    the_var = prefix+'pf_avg'
    store_data, the_var, orbit_times, data
    add_setting, the_var, /smart, dictionary($
        'display_type', 'scalar', $
        'unit', unit, $
        'short_name', 'RB-'+strupcase(probe))
endforeach


time_step = constant('secofday')
common_times = make_bins(time_range, time_step)
apogee_mlt = (get_var_data('rbspa_apogee_mlt', at=common_times)+$
    get_var_data('rbspb_apogee_mlt', at=common_times))*0.5
mlt_range = [-1,1]*7.5
index = where_pro(apogee_mlt, '[]', mlt_range)
mlt_time_ranges = time_to_range(common_times[index], time_step=time_step)
durations = (mlt_time_ranges[*,1]-mlt_time_ranges[*,0])
index = where(durations ge 4*30*time_step)
mlt_time_ranges = mlt_time_ranges[index,*]
dmlt = 3
mlt_constants = make_bins(mlt_range, dmlt*2, /inner)

the_var = 'apogee_mlt'
store_data, the_var, common_times, [$
    [get_var_data('rbspa_apogee_mlt', at=common_times)], $
    [get_var_data('rbspb_apogee_mlt', at=common_times)]]
add_setting, the_var, /smart, dictionary($
    'display_type', 'list', $
    'unit', 'hr')
options, the_var, 'colors', sgcolor(['red','blue'])
options, the_var, 'labels', 'RB-'+strupcase(probes)
options, the_var, 'constant', mlt_constants
options, the_var, 'ytickv', mlt_constants
options, the_var, 'yminor', dmlt
options, the_var, 'yticks', n_elements(mlt_constants)-1
options, the_var, 'yrange', mlt_range
options, the_var, 'ystyle', 1


the_var = 'dst'
pflux_survey_load_data, the_var, probe='b'
rename_var, 'rbspb_dst', to=the_var
options, the_var, 'constant', [0,-40]
options, the_var, 'yrange', [-150,50]
options, the_var, 'ytickv', [0,-100]
options, the_var, 'yticks', 1
options, the_var, 'yminor', 5

vars = tnames('rbsp?_pf_avg')
options, vars, 'ylog', 1
options, vars, 'yrange', [0.2,20]
get_data, 'rbspa_pf_avg', times, pfa
pfb = get_var_data('rbspb_pf_avg', at=times)
the_var = 'pf_avg'
store_data, the_var, times, [[pfa],[pfb]]
add_setting, the_var, /smart, dictionary($
    'display_type', 'list', $
    'unit', 'mW/m!U2!N')
options, the_var, 'colors', sgcolor(['red','blue'])
options, the_var, 'labels', 'RB-'+strupcase(probes)
options, the_var, 'ylog', 1
options, the_var, 'yrange', [0.2,20]

nmlt_range = 1

project = pflux_survey_load_project()
plot_file = join_path([project.plot_dir,'fig_pflux_vs_dst.pdf'])
if keyword_set(test) then plot_file = test
sgopen, plot_file, xsize=8, ysize=4
margins = [12,4,5,2]
all_poss = sgcalcpos(nmlt_range, ypad=4, margins=margins, xchsz=xchsz, ychsz=ychsz)
xticklen_chsz = -0.15
yticklen_chsz = -0.30

npanel = 3
for range_id=0, nmlt_range-1 do begin
    poss = sgcalcpos(npanel, position=all_poss[*,range_id], ypans=[2,2,1])
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    vars = ['dst','pf_avg','apogee_mlt']
    options, vars, 'xticklen', xticklen
    options, vars, 'yticklen', yticklen
    
    the_time_range = mlt_time_ranges[range_id,*]
    tplot, vars, position=poss, /noerase, trange=the_time_range
    
    labels = letters(npanel)+'. '+['Dst','Orb.Avg.S','Apogee!C    MLT']
    for ii=0,npanel-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,/normal, labels[ii]
    endfor
endfor
if keyword_set(test) then stop
sgclose

end
