
pro _2013_0607_load_burst_data3, filename=filename, $
    reload_burst_field=reload_burst_field, $
    reload_vb1=reload_vb1

    
    id = '2013_0607'
    tvar = id+'_event_info'
    if tnames(tvar) eq '' then _2013_0607_load_data3
    get_data, tvar, tmp, einfo
    time_range = einfo.time_range_short
    probes = einfo.probes
    pres = einfo.pres
    
    uvw = ['u','v','w']
    rgb = sgcolor(['red','green','blue'])
    
    update_datfn = 0
    

;---Check inputs.
    ; check datfn.
    load = 0
    if n_elements(filename) eq 0 then load = 1 else if file_test(filename) eq 0 then load = 1
    if load eq 1 then begin
        rootdir = sdiskdir('GoogleDrive')+'/My Drive/works/works/rbsp_psbl1'
        datdir = rootdir+'/data'
        if file_test(datdir,/directory) eq 0 then file_mkdir, datdir
        datfn = datdir+'/2013_0607_burst_data_v3.tplot'
    endif else datfn = filename
    if file_test(datfn) eq 1 then tplot_restore, filename=datfn
    
    
;---Load burst data.
    burst_data_rates = 1d/[4096d,1024]  ; in sec.
    
    load = keyword_set(reload_burst_field)
    tvar = ['eb1_uvw','bb1_uvw']
    field_vars = []
    foreach pre0, pres do field_vars = [field_vars,pre0+tvar]
    foreach tvar, field_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        foreach probe, probes, ii do begin
            rbsp_read_b1_field, time_range, probe=probe, prefix=pres[ii], datarate=burst_data_rates[ii]
        endforeach
        update_datfn = 1
    endif
    
;---Load probe potentials.
    nboom = einfo.nboom
    boom_labels = einfo.boom_labels
    boom_colors = einfo.boom_colors
    load = keyword_set(reload_vb1)
    tvar = ['vb1']
    vb1_vars = []
    foreach pre0, pres do vb1_vars = [vb1_vars,pre0+tvar]
    foreach tvar, vb1_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        foreach probe, probes, ii do begin
            rbsp_load_efw_waveform, probe=probe, type='calibrated', $
                datatype='vb1', /downloadonly, files=fn, ignore_filedate=1
            
            ; read data only within the wanted time.
            cdfid = cdf_open(fn)
            cdf_control, cdfid, variable = 'epoch', get_var_info = vinfo, /zvariable
            nrec = vinfo.maxrec+1
            cdf_varget, cdfid, 'epoch', et16, rec_start = 0, rec_count = nrec
            et16 = transpose(et16)
            ; convert epoch16 to ut sec.
            uts = real_part(et16)+imaginary(et16)*1d-12 - 62167219200D
            et16 = 0    ; free memory.
            ; find the record range.
            idx = where(uts ge time_range[0] and uts le time_range[1], nrec)
            rec0 = idx[0]
            ; read variables, epoch, vb1.
            epname = 'epoch'
            vname = 'vb1'
            cdf_varget, cdfid, epname, uts, rec_start = rec0, rec_count = nrec & uts = transpose(uts) & uts = real_part(uts)+imaginary(uts)*1d-12 - 62167219200D
            cdf_varget, cdfid, vname, vb1, rec_start = rec0, rec_count = nrec & vb1 = transpose(vb1)
            cdf_close, cdfid

            ; get calibration parameters.
            cp0 = rbsp_efw_get_cal_params(time_range[0])
            case strlowcase(probe) of
                'a': cp = cp0.a
                'b': cp = cp0.b
            endcase
            gain = cp.adc_gain_vdc
            offset = cp.adc_offset_vdc
            vb1 = double(vb1)
            for i = 0,nboom-1 do vb1[*,i] = (vb1[*,i]-offset[i])*gain[i]
            store_data, pres[ii]+'vb1', uts, vb1, limits={$
                ytitle:'(V)', labels:'V'+boom_labels, colors:boom_colors}
            
            ; interpolate to uniform times.
            tvar = ['vb1']
            foreach var, tvar do uniform_time, pres[ii]+tvar, burst_data_rates[ii]
        endforeach
        update_datfn = 1
    endif


;---Save changes.
    if update_datfn then begin
        vars_to_save = [field_vars, vb1_vars]
        lprmsg, 'Saving data to '+datfn+' ...'
        tplot_save, vars_to_save, filename=datfn
    endif

end


_2013_0607_load_burst_data3, reload_vb1=1
end