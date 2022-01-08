;+
; Read the time range of sections shifted in time tag, and the correction.
; Data adopted from the txt file generated by rbsp_efw_read_l1_time_tag_correction_gen_file.
; 
; The tplot_var is 'rbspx_l1_time_tag_correction', to use
;    get_data, 'rbspa_l1_time_tag_correction', start_times, time_ranges, corrections
;    nsection = n_elements(corrections)
;   
;    foo_l1_efw_data = 'rbspa_efw_esvy'
;    get_data, foo_l1_efw_data, times, data
;    for ii=0, nsection-1 do begin
;        tmp = where(times ge time_ranges[ii,0] and times le time_ranges[ii,1], count)
;        if count eq 0 then continue
;        ; Have to find the closest time, otherwise the index can be 1 record off.
;        if min(times) ge correction_time_ranges[ii,0] then i0 = 0 else begin
;            index = min(times-correction_time_ranges[ii,0], /absolute, i0)
;        endelse
;        if max(times) le correction_time_ranges[ii,1] then i1 = n_elements(times) else begin
;            index = min(times-correction_time_ranges[ii,1], /absolute, i1)
;        endelse
;        times[i0:i1-1] += corrections[ii]
;    endfor
;    store_data, foo_l1_efw_data, times, data
;-

pro rbsp_efw_read_l1_time_tag_correction, probe=probe

    info = [$
        'a','2014-01-05/12:42:34.450798','2014-01-05/13:27:05.419235','-1.0000019', $
        'a','2014-04-23/18:17:28.683166','2014-04-23/18:57:11.651718','-1.0000010', $
        'a','2014-04-24/15:26:31.677452','2014-04-24/16:06:46.645874','-0.9999862', $
        'a','2014-05-22/06:31:51.498458','2014-05-22/07:10:14.466918','-0.9999862', $
        'a','2014-06-16/00:00:07.590774','2014-06-16/04:52:24.342498',' 1.0000150', $
        'a','2014-06-20/00:00:02.727035','2014-06-20/08:09:23.317199',' 0.9999990', $
        'a','2014-07-08/07:28:46.208732','2014-07-08/08:05:33.177345','-1.0000012', $
        'a','2014-07-09/03:09:01.203872','2014-07-09/03:46:20.172317','-0.9999859', $
        'a','2014-07-09/22:55:24.198844','2014-07-09/23:32:43.167533','-1.0000148', $
        'a','2014-07-11/14:26:34.189186','2014-07-11/15:03:05.157585','-1.0000000', $
        'a','2016-08-09/23:59:51.192207','2016-08-10/00:36:22.160964','-1.0000000', $
        'a','2016-10-19/08:28:43.235839','2016-10-19/09:02:34.204597','-1.0000000', $
        'a','2016-10-20/01:57:14.236579','2016-10-20/02:31:05.205253','-1.0000000', $
        'a','2016-10-20/19:27:37.236976','2016-10-20/20:00:56.205871','-1.0000000', $
        'a','2016-10-21/12:57:28.237754','2016-10-21/13:30:15.206344','-1.0000150', $
        'a','2016-10-22/06:26:47.238235','2016-10-22/06:58:46.207061','-1.0000160', $
        'a','2016-10-22/23:55:50.238731','2016-10-23/00:00:05.234603','-0.9999809', $
        'a','2016-10-29/13:22:21.244438','2016-10-29/13:55:24.213134','-1.0000000', $
        'a','2016-11-02/05:12:40.247650','2016-11-02/05:47:03.216316','-1.0000010', $
        'a','2016-11-02/22:52:07.248321','2016-11-02/23:27:34.217048','-0.9999850', $
        'a','2016-11-11/07:29:00.256393','2016-11-11/08:04:11.225013','-0.9999840', $
        'a','2016-11-14/22:25:12.260253','2016-11-14/23:13:43.228935','-1.0000160', $
        'a','2016-11-20/00:00:04.239898','2016-11-20/00:09:07.235038','-1.0000141', $
        'a','2016-11-21/06:48:35.267776','2016-11-21/07:49:22.236633','-1.0000000', $
        'a','2016-11-24/02:32:49.271362','2016-11-24/03:41:36.240249','-1.0000000', $
        'a','2019-06-06/00:00:29.803977','2019-06-06/00:03:40.801582','-1.0000138', $
        'b','2014-04-22/19:58:55.036964','2014-04-22/22:49:02.060295','-1.0000050', $
        'b','2014-06-15/00:00:21.690292','2014-06-18/00:00:23.089233',' 0.9999940', $
        'b','2016-04-07/23:58:59.704216','2016-04-08/02:10:10.728927','-0.9999890']

    ndim = 4
    nrec = n_elements(info)/ndim
    infos = transpose(reform(info, ndim,nrec))

    index = where(infos[*,0] eq probe, nrec)
    infos = infos[index,*]

    tformat = 'YYYY-MM-DD/hh:mm:ss.ffffff'
    time_ranges = time_double(infos[*,1:2], tformat=tformat)

    corrections = double(infos[*,3])

    prefix = 'rbsp'+probe+'_'
    var = prefix+'l1_time_tag_correction'
    store_data, var, time_ranges[*,0], time_ranges, corrections

end

rbsp_efw_read_l1_time_tag_correction, probe='a'
rbsp_efw_read_l1_time_tag_correction, probe='b'
end