
pro psbl_de_trim_b1_data, utr, probe, type = dattype

    if n_elements(utr) ne 2 then message, 'no time range ...'
    if n_elements(probe) ne 1 then message, 'no probe ...'
    if n_elements(dattype) eq 0 then dattype = 'vb1'
    
    timespan, utr[0], utr[1]-utr[0], /second
    
    pre0 = 'rbsp'+probe+'_'
    
    ; download b1 data file.
    rbsp_load_efw_waveform, probe = probe, type = 'calibrated', $
        datatype = dattype, /downloadonly, files = fn, ignore_filedate = 1
    
    ; read data only within the wanted time.
    cdfid = cdf_open(fn)
    cdf_control, cdfid, variable = 'epoch', get_var_info = vinfo, /zvariable
    nrec = vinfo.maxrec+1
    cdf_varget, cdfid, 'epoch', et16, rec_start = 0, rec_count = nrec
    et16 = transpose(et16)
    ; convert epoch16 to ut sec.
    uts = real_part(et16)+imaginary(et16)*1d-12 - 62167219200D
    et16 = 0    ; free memory.
    
    idx = where(uts ge utr[0] and uts le utr[1], nrec)
    if nrec eq 0 then begin
        cdf_close, cdfid
        return
    endif
    rec0 = idx[0]
    
    ; read variables, epoch, vb1.
    case dattype of
        'vb1': begin
            epname = 'epoch'
            vname = 'vb1'
            end
        'mscb': begin
            epname = 'epoch'
            vname = 'mscb'
            end
    endcase
    cdf_varget, cdfid, epname, uts, rec_start = rec0, rec_count = nrec & uts = transpose(uts) & uts = real_part(uts)+imaginary(uts)*1d-12 - 62167219200D
    cdf_varget, cdfid, vname, vb1, rec_start = rec0, rec_count = nrec & vb1 = transpose(vb1)

    cdf_close, cdfid
    
    
    ; get calibration parameters.
    cp0 = rbsp_efw_get_cal_params(utr[0])
    case strlowcase(probe) of
        'a': cp = cp0.a
        'b': cp = cp0.b
        else: dprint, 'Invalid probe name. Calibration aborted.'
    endcase
    gain = cp.adc_gain_vdc
    offset = cp.adc_offset_vdc
    vb1 = double(vb1)
    for i = 0, 5 do vb1[*,i] = (vb1[*,i]-offset[i])*gain[i]
stop    
    
    store_data, pre0+'vb1', uts, vb1, limits = $
        {ytitle:'V B1!C(V)', labels:'V'+['1','2','3','4','5','6'], colors:[1,2,3,4,5,6]}
    uts = 0
    vb1 = 0
    

    ; download mscb
    rbsp_load_efw_waveform, probe = probe, type = 'calibrated', $
        datatype = ['mscb'], /downloadonly, files = fn, ignore_filedate = 1
stop
    ; read data only within the wanted time.
    cdfid = cdf_open(fn)
    cdf_control, cdfid, variable = 'epoch', get_var_info = vinfo, /zvariable
    nrec = vinfo.maxrec+1
    cdf_varget, cdfid, 'epoch', et16, rec_start = 0, rec_count = nrec
    et16 = transpose(et16)
    ; convert epoch16 to ut sec.
    uts = real_part(et16)+imaginary(et16)*1d-12 - 62167219200D
    et16 = 0    ; free memory.
    
    idx = where(uts ge utr[0] and uts le utr[1], nrec)
    if nrec eq 0 then begin
        cdf_close, cdfid
        return
    endif
    rec0 = idx[0]
    
    ; read variables, epoch, vb1.
    cdf_varget, cdfid, 'epoch', uts, rec_start = rec0, rec_count = nrec & uts = transpose(uts) & uts = real_part(uts)+imaginary(uts)*1d-12 - 62167219200D
    cdf_varget, cdfid, '', vb1, rec_start = rec0, rec_count = nrec & vb1 = transpose(vb1)
    
    cdf_close, cdfid
    
    
end


utr = time_double(['2013-04-14/08:02','2013-04-14/08:06'])
utr = time_double(['2013-06-07/04:53','2013-06-07/04:58'])
probes = ['a','b']
foreach tprobe, probes do psbl_de_trim_b1_data, utr, tprobe
end