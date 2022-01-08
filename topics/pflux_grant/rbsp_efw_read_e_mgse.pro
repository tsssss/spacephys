;+
; Preprocess and load e_mgse to memory.
; 
; keep_spin_axis=. Set to keep spin axis data. Set to it to 0 by default.
;-

pro rbsp_efw_read_e_mgse, time_range, probe=probe, keep_spin_axis=keep_spin_axis

    prefix = 'rbsp'+probe+'_'
    the_var = prefix+'e_mgse'
;    if ~check_if_update(the_var, time_range, dtime=60) then return


;---New codes, as of 2021/07/05.
    e_uvw_var = prefix+'e_uvw'
    rbsp_efw_phasef_read_e_uvw, time_range, probe=probe, keep_spin_axis=keep_spin_axis
    get_data, e_uvw_var, times
    if n_elements(times) lt 10 then begin
        errmsg = 'No E UVW ...'
        return
    endif
    
;    common_time_step = 1d/16
;    common_times = make_bins(time_range, common_time_step)
;    interp_time, e_uvw_var, common_times
    
    e_uvw = get_var_data(e_uvw_var, times=common_times)
    rbsp_read_q_uvw2gse, time_range, probe=probe
    if ~keyword_set(keep_spin_axis) then e_uvw[*,2] = 0
    e_mgse = cotran(e_uvw, common_times, probe=probe, 'uvw2mgse')
    store_data, the_var, common_times, e_mgse
    add_setting, the_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'mV/m', $
        'short_name', 'E', $
        'coord', 'MGSE', $
        'coord_labels', constant('xyz') )

end


time_range = time_double(['2013-02-22','2013-02-24'])   ; Data gap between days.
time_range = time_double(['2016-10-27','2016-10-28'])   ; Data gap within a day.
time_range = time_double(['2014-04-08','2014-04-10'])   ; Data gap within a day.
time_range = time_double(['2015-09-30','2015-10-01'])   ; Jump in E MGSE, caused by data gap.
time_range = time_double(['2016-10-21','2016-10-22'])   ; Jump in E MGSE.
time_range = time_double(['2016-10-22','2016-10-23'])   ; Jump in E MGSE.
probe = 'a'
;time_range = time_double(['2012-10-06','2012-10-07'])
;probe = 'a'
rbsp_efw_read_e_mgse, time_range, probe=probe
end
