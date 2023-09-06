;+
; For theta_var and a given time_range, calculate:
;   1. theta_smooth. A smoothed version of theta.
;   2. theta_stddev. The standard deviation of theta.
;   3. dtheta. A smoothed dtheta.
;   4. dtheta_stddev. The standard deviation of dtheta.
; Save these data to theta_var+'_smooth_combo'
;
; theta_var. A string for the tplot var.
; time_range. 2 doubles for the time range.
; smooth_window=. A number in sec for the smooth window.
; stat_ratio=. A number in [0,1] to set the range of points to do stats.
;-

pro azim_df_smooth_theta, theta_var, time_range, $
    smooth_window=smooth_window, stat_ratio=stat_ratio

    if n_elements(smooth_window) eq 0 then smooth_window = 120. ; sec.
    if n_elements(stat_ratio) eq 0 then stat_ratio = 0.8

    theta = get_var_data(theta_var, in=time_range, times=times)

    ; Prepare low-res times.
    boxcar_boundary_times = smkarthm(time_range[0], time_range[1], smooth_window, 'dx')
    nboxcar = n_elements(boxcar_boundary_times)-1
    boxcar_center_times = boxcar_boundary_times[0:nboxcar-1]+smooth_window*0.5

    ; Calculate theta_smooth, theta_stddev.
    theta_smooth = fltarr(nboxcar)
    theta_stddev = fltarr(nboxcar)
    for ii=0, nboxcar-1 do begin
        index = where_pro(times, '[]', boxcar_boundary_times[ii:ii+1], count=count)
        if count eq 0 then continue
        the_theta = theta[index]
        ; Remove nans.
        index = where(finite(the_theta), count)
        if count lt 2 then continue
        the_theta = the_theta[index]
        ; Use points around mean to reduce jitter.
        index = sort(abs(the_theta-mean(the_theta)))
        index = index[0:count*stat_ratio]
        the_theta = the_theta[index]
        ; Do statistics.
        theta_smooth[ii] = mean(the_theta)
        theta_stddev[ii] = stddev(the_theta)
    endfor
    dtheta = deriv(boxcar_center_times, theta_smooth)

    theta_smooth = interpol(theta_smooth, boxcar_center_times, times, /quadratic)
    theta_stddev = interpol(theta_stddev, boxcar_center_times, times, /quadratic)
    dtheta = interpol(dtheta, boxcar_center_times, times, /quadratic)

    ; Calculate dtheta and its stddev.
    dtheta_stddev = fltarr(nboxcar)+!values.f_nan
    for ii=0, nboxcar-1 do begin
        index = where_pro(times, '[]', boxcar_boundary_times[ii:ii+1], count=count)
        if count lt 2 then continue
        the_dtheta = dtheta[index]
        dtheta_stddev[ii] = stddev(the_dtheta, /nan)
    endfor
    min_dtheta = stddev(dtheta_stddev, /nan)
    index = where(finite(dtheta_stddev,/nan), count)
    if count ne 0 then dtheta_stddev[index] = min_dtheta
    dtheta_stddev = smooth(dtheta_stddev>min_dtheta, 7, /edge_truncate)
    dtheta_stddev = interpol(dtheta_stddev, boxcar_center_times, times)

;---Save data.
    theta_smooth_combo = theta_var+'_smooth_combo'
    store_data, theta_smooth_combo, times, [[theta_smooth],[theta_stddev],[dtheta],[dtheta_stddev]]

;    theta_combo_var = theta_var+'_combo'
;    store_data, theta_combo_var, times, [[theta],[theta_smooth],[theta_stddev],[-theta_stddev]], $
;        limits={ytitle:'(deg)',labels:tex2str('theta')+' '+['orig','smooth','stddev',''], colors:sgcolor(['silver','black','tan','tan']), constant:0}
;    dtheta_combo_var = theta_var+'_dtheta_combo'

end
