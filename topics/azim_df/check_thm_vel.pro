;+
; Check themis esa ion velocity.
;-

time_range = time_double(['2014-08-28/10:00','2014-08-28/10:30'])
df_list = list()
foreach probe, probes do begin
    dfs = azim_df_read_df(time_range, probe='th'+probe)
    df_list.add, dfs, /extract
endforeach

probes = ['a','d','e']
;foreach probe, probes do themis_read_ion_vel, time_range, probe=probe
tplot, 'th?_u_gsm', trange=time_range

v_list = list()
foreach df, df_list do begin
    prefix = df.probe+'_'
    v_gsm = get_var_data(prefix+'u_gsm', in=df.obs_time+[-1,1]*5*60, times=times)
    v_gsm = total(v_gsm,1)/n_elements(times)
    v_sm = cotran(v_gsm, df.obs_time, 'gsm2sm')
    v_list.add, v_sm
    timebar, df.obs_time
endforeach

foreach v_sm, v_list do print, v_sm
end