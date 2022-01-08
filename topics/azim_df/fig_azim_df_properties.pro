;+
; Plot properties of the DFs in the azimuthal DF events.
;-

events = azim_df_load_azim_events(project=project)


;keys = events[0].df_list[0].keys()
;property_info = dictionary()
;foreach key, keys do begin
;    the_list = list()
;    foreach event, events do begin
;        foreach df, event.df_list do the_list.add, df[key]
;    endforeach
;    property_info[key] = the_list.toarray()
;endforeach
;
;mlts = azim_df_calc_pseudo_mlt(property_info.arrival_r_sm)
;rxys = snorm(property_info.arrival_r_sm[*,0:1])
;heights = property_info.height
;widths = property_info.width

event_times = list()
foreach event, events do begin
    event_times.add, mean(event.time_range) 
endforeach
event_times = event_times.toarray()

event_durations = list()
foreach event, events do begin
    data = list()
    foreach df, event.df_list do data.add, df.arrival_time
    df_times = data.toarray()
    event_durations.add, total(minmax(df_times)*[-1,1])
endforeach
event_durations = event_durations.toarray()

rxys = list()
foreach event, events do begin
    data = list()
    foreach df, event.df_list do data.add, df.arrival_r_sm
    r_sms = data.toarray()
    rxys.add, minmax(snorm(r_sms[*,0:1]))
endforeach
rxys = rxys.toarray()

mlts = list()
foreach event, events do begin
    data = list()
    foreach df, event.df_list do data.add, df.arrival_r_sm
    r_sms = data.toarray()
    mlt = azim_df_calc_pseudo_mlt(r_sms)
    mlts.add, sign(mlt[0])*minmax(abs(mlt))
endforeach
mlts = mlts.toarray()

;plot, [0,30], [10,1000], /nodata, ylog=1
vmags = list()
omegas = list()
timing_key = project.timing_key
km_s2deg_min = 1d/constant('re')*60*constant('deg')
foreach event, events do begin
    timing_info = event.timing_info[timing_key]
    combos = timing_info.triad_timing
    good_triad_flags = timing_info.good_triad_flags
    vmag = list()
    omega = list()
    rxy = list()
    foreach key, combos.keys(), ii do begin
        combo = combos[key]
        if good_triad_flags[ii] eq 0 then continue
        the_rxy = snorm(combo.center_rsm[0:1])
        rxy.add, the_rxy
        vmag.add, combo.v_mag
        omega.add, combo.v_mag/the_rxy*km_s2deg_min
    endforeach
    vmag = vmag.toarray()
    rxy = rxy.toarray()
;    plots, rxy, vmag, psym=-1
;    stop
    vmags.add, vmag
    omegas.add, omega.toarray()
endforeach

vmag_mean = fltarr(events.length)
vmag_stddev = fltarr(events.length)
foreach data, vmags, ii do begin
    vmag_mean[ii] = mean(data)
    vmag_stddev[ii] = stddev(data)
endforeach

omega_mean = fltarr(events.length)
omega_stddev = fltarr(events.length)
foreach data, omegas, ii do begin
    omega_mean[ii] = mean(data)
    omega_stddev[ii] = stddev(data)
endforeach




omega_fit = list()
r2 = list()
foreach event, events do begin
    timing_info = event.timing_info[timing_key]
    r2.add, timing_info.linear_fit_info.r_square
    omega_fit.add, timing_info.linear_fit_info.omega_fit
endforeach
omega_fit = omega_fit.toarray()
r2 = r2.toarray()


omega_range = [0,10]
xrange = omega_range
yrange = omega_range
plot, xrange, yrange, /nodata, /iso
foreach event, events, ii do begin
    txs = omega_mean[ii]+[-1,1]*omega_stddev[ii]*0.5
    tys = abs(omega_fit[ii])+[0,0]
    plots, txs, tys, /data
    plots, mean(txs), tys[0], psym=8
endforeach
plots, xrange, yrange, linestyle=1
stop

plot, [0,30], [0,20], /nodata
widths = list()
re = constant('re')
foreach event, events, event_id do begin
    width = list()
    rxy = list()
    foreach df, event.df_list do begin
        width.add, df.width
        rxy.add, snorm(df.arrival_r_sm[0:1])
    endforeach
    rxy = rxy.toarray()
    width = width.toarray()
    index = sort(rxy)
    rxy = rxy[index]
    width = width[index]
;    plots, rxy, width, psym=-1
;    stop
    widths.add, width
endforeach
width_in_re = fltarr(events.length)
width_in_deg = fltarr(events.length)
foreach width, widths, ii do begin
    the_width = median(width)
    width_in_re[ii] = the_width*vmag_mean[ii]/re
    width_in_deg[ii] = the_width*omega_mean[ii]/60
endforeach


ae = list()
foreach event, events do begin
    ae.add, event.max_ae
endforeach
ae = ae.toarray()


stop

end
