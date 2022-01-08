;+
; Preprocess B field for the whole mission, to:
;   1. remove wobble at harmonics of spin freq.
;   2. remove spikes.
;   3. change fillval -99999 to data gap.
;-

pro pflux_grant_read_b_gsm_gen_file, time, probe=probe, filename=file;, pad_time=pad_time

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    valid_range = time_double(['2012-09-29','2015-10-03'])
    if n_elements(file) eq 0 then begin
        errmsg = handle_error('No output file ...')
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif

    if n_elements(time) eq 0 then begin
        errmsg = handle_error('No input time ...')
        return
    endif


;---Settings.
    secofday = constant('secofday')
    date = time[0]-(time[0] mod secofday)
    day_time_range = date+[0,secofday]
    pad_time = 300. ; pad to remove boundary effect due to wavelet transform.
    time_range = day_time_range+[-1,1]*pad_time
    project = pflux_grant_load_project()
    common_time_step = project.common_time_step
    common_times = make_bins(time_range, common_time_step)
    ncommon_time = n_elements(common_times)
    ndim = 3
    prefix = 'rbsp'+probe+'_'
    xyz = constant('xyz')
    fillval = !values.f_nan
    min_ntime = 1800d*64    ; 0.5 hour of data.

;---Load B field.
    b_gsm_var = prefix+'b_gsm'
    b_uvw_var = prefix+'b_uvw'
    rbsp_read_emfisis, time_range, probe=probe, id='l2%magnetometer'

    get_data, b_uvw_var, bfield_times
    nbfield_time = n_elements(bfield_times)
    if nbfield_time le min_ntime then begin
        b_gsm = fltarr(ncommon_time,ndim)+fillval
        bfield_times = common_times
    endif else begin
        ; Fix wobble.
        pflux_grant_fix_b_uvw, time_range, probe=probe

        ; Convert to GSM.
        rbsp_read_q_uvw2gse, time_range, probe=probe
        b_uvw = get_var_data(b_uvw_var, times=bfield_times)
        b_gsm = cotran(b_uvw, bfield_times, 'uvw2gsm', probe=probe)
    endelse

    store_data, b_gsm_var, bfield_times, b_gsm
    add_setting, b_gsm_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz )
    if n_elements(bfield_times) ne ncommon_time then interp_time, b_gsm_var, common_times




;---Save data to file.
    index = lazy_where(common_times, '[)', day_time_range)
    times = common_times[index]
    b_gsm = get_var_data(b_gsm_var)
    b_gsm = float(b_gsm[index,*])


    cdf_id = (file_test(file))? cdf_open(file): cdf_create(file)

    time_var = 'ut_sec'
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var


    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz)
    cdf_save_var, b_gsm_var, value=b_gsm, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=b_gsm_var

    cdf_close, cdf_id

end

probe = 'b'
;time = time_double('2012-12-04')    ; gap.
;time = time_double('2012-10-02')    ; spikes.
time = time_double('2014-08-28')    ; spin tone.
time = time_double('2012-10-01')    ; storm.

probe = 'a'
time = time_double('2012-10-13')    ; storm.
time = time_double('2012-10-08')    ; storm.
time = time_double('2013-01-01')    ; Drift.
time = time_double('2012-10-11')    ; spikes.

file = join_path([homedir(),'test.cdf'])
pflux_grant_read_b_gsm_gen_file, time, probe=probe, filename=file

stop

pad_times = [300,3600]
npad_time = n_elements(pad_times)
colors = constant('rgb')
labels = string(pad_times,format='(I0)')
ndim = 3
xyz = constant('xyz')
run_times = fltarr(npad_time)
prefix = 'rbsp'+probe+'_'

foreach pad_time, pad_times, ii do begin
    file = join_path([homedir(),'test_'+string(pad_time,format='(I0)')+'_sec.cdf'])
    if file_test(file) then file_delete, file
    del_data, '*'
    tic
    pflux_grant_read_b_gsm_gen_file, time, probe=probe, filename=file, pad_time=pad_time
    run_times[ii] = toc()
endforeach



foreach pad_time, pad_times, ii do begin
    file = join_path([homedir(),'test_'+string(pad_time,format='(I0)')+'_sec.cdf'])
    var = prefix+'b_gsm'
    cdf_load_var, var, filename=file
    get_data, var, times, b_gsm
    if ii eq 0 then begin
        ntime = n_elements(times)
        b = fltarr(ntime,npad_time,ndim)
        b[*,ii,*] = b_gsm
    endif else begin
        b[*,ii,*] = b_gsm
    endelse
endforeach

for ii=0,ndim-1 do begin
    store_data, prefix+'b'+xyz[ii], times, b[*,*,ii], limits={colors:colors, labels:labels}
endfor
store_data, prefix+'db', times, b[*,1,2]-b[*,0,2]
sgopen
tplot, 'rbspb_'+['b'+xyz,'db'], trange=time+[0,86400]
print, run_times

end
