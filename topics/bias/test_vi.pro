


pro test_vi, tr, rbx, probe, reload = reload, vrange = vrng

    ; parameters.
;    rbx = 'a'
;    probe = 3
;    date = '2014-12-17'
;    date = '2013=09-11'
    
    dt = 86400d
 

    ; s/c, can be 'a' or 'b'.
;    rbx = 'a'
    rbspx = 'rbsp'+rbx
    prefix = rbspx+'_'
    
    ; probes, can be [1,2], [3,4], [5,6].
;    if n_elements(probe) eq 0 then probe = [1,2,3,4,5,6]
    suffix = string(probe,format='(I0)')
    pair = string(22*(probe/2)+12,format='(I2)')
    
    ; time range to get V-I curve.
    timespan, tr[0], tr[1]-tr[0], /second
    
    
    ; load data or not.
    load = 0
    
    ; no data have been loaded.
    vars = prefix+['V'+suffix,['guard','usher','ibias']+suffix, $
        'efw_vsvy','efw_esvy']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], data = data
        ; not a combo or structure.
        if n_elements(data) eq 1 and size(data,/type) ne 8 then load = 1
    endfor
    
    ; new time is set.
    get_data, vars[0], tmp, data
    if n_elements(tmp) ne 1 then begin
        if min(tmp) gt max(tr) or max(tmp) lt min(tr) then load = 1
    endif
    
    ; force to reload.
    if keyword_set(reload) then load = 1
    
    if load eq 1 then begin
        rbsp_efw_init
        
        ; load bias current.
        ; rbspx_IEFI_[GUARD,IBIAS,USHER][123456].
        type = 'hsk_beb_analog'
        files = rbspx+'/l1/'+type+'/YYYY/'+rbspx+'_l1_'+type+'_YYYYMMDD_v*.cdf'
        files = file_dailynames(file_format=files,trange=tr)
        files = file_retrieve(files, /last_version, _extra=!rbsp_efw)
        cdf2tplot, file = files, /get_support_data, /convert_int1_to_int2, prefix = prefix
    
        vars = prefix+['CONFIG0','IBEB_TEMP','SPARE','gdh_*','ga_*','ccsds_*','packet_*']
        store_data, vars, /delete
        
        ; convert to step style.
        for i = 0, 3-1 do begin
            tsufx = [string(i*2+1,format='(I0)'),string(i*2+2,format='(I0)')]
            tcomp = strjoin(tsufx,'')
            
            ; find step time, convert to step style.
            vars = tnames(prefix+'IEFI_*'+tsufx)
            for j = 0, n_elements(vars)-1 do begin
;                bias_param_info, vars[j], trange = tr, datarate = dr0
                tvar = strlowcase(prefix+strmid(vars[j],strpos(vars[j],'_',/reverse_search)+1))
                bias_param_to_step, vars[j], newname = tvar
            endfor
        endfor

        ; load probe potential.
        ; rbspx_efw_vsvy.
        rbsp_load_efw_waveform, probe = rbx, type = 'calibrated', datatype = ['esvy','vsvy']
        vars = prefix+['efw_esvy_*','efw_vsvy_*']
        store_data, vars, /delete

        ; split probe potentials into components.
        vars = prefix+'efw_vsvy'
        get_data, vars, t0, vsvy
        for i = 0, 6-1 do begin
            tvar = prefix+'V'+string(i+1,format='(I0)')
            store_data, tvar, t0, vsvy[*,i]
            options, tvar, 'ysubtitle', '[V]'
            options, tvar, 'colors', i
            options, tvar, 'labels', 'V'+string(i+1,format='(I0)')
        endfor
        
        for i = 0, 3-1 do begin
            tvar = prefix+'Vsc'+string(22*i+12,format='(I2)')
            store_data, tvar, t0, 0.5*(vsvy[*,2*i]+vsvy[*,2*i+1])
            options, tvar, 'ysubtitle', '[V]'
            options, tvar, 'colors', i*2+2
            options, tvar, 'labels', 'Vsc'+string(22*i+12,format='(I2)')
        endfor
    endif
    
    ; 'rbspx_IEFI_[IBIAS,GAURD,USHER][123456]' save the original parameters.
    ; 'rbspx_[ibias,guard,usher][123456]' save the stepped version.


    ; read voltage-current relation.    
    tvar = prefix+'Vsc'+pair
    get_data, tvar, t1, vscs
    
    tvar = prefix+'ibias'+suffix
    get_data, tvar, t0, ibias
    tvar = prefix+'usher'+suffix
    get_data, tvar, t0, usher
    tvar = prefix+'guard'+suffix
    get_data, tvar, t0, guard
    
    idx = where(t0 ge tr[0] and t0 le tr[1])
    t0 = t0[idx]
    ibias = ibias[idx]
    usher = usher[idx]
    guard = guard[idx]
    
    ; find the bias period.
    biasutr = [0d,0]
    dibias = 5  ; nA.
    dt = 40     ; sec.
    minibias = min(ibias)
    idx = where(ibias le minibias+dibias)
    maxibias = max(ibias[min(idx):max(idx)])
    biasutr[0] = t0[idx[0]]-dt
    idx = where(ibias ge maxibias-dibias)
    biasutr[1] = t0[idx[-1]]+dt
    
    idx = where(t0 ge biasutr[0] and t0 le biasutr[1])
    t0 = t0[idx]
    ibias = ibias[idx]
    usher = usher[idx]
    guard = guard[idx]


    usher0 = smode(usher)
    guard0 = smode(guard)
    err = 0.5   ; V.
    
    vars = prefix+['Vsc'+pair,'ibias'+suffix]
    tplot, vars, trange = biasutr
        
    nrec = n_elements(t0)
    vs = dblarr(nrec)
    for i = 0, nrec-2 do begin
;        if abs(usher[i]-usher0) gt err then continue
;        if abs(guard[i]-guard0) gt err then continue
        idx = where(t1 ge t0[i] and t1 lt t0[i+1])
        vs[i] = mean(vscs[idx])
    endfor
    
    
;    idx = where(vs ne 0)
;    vs = vs[idx]
;    ibias = ibias[idx]
;    usher = usher[idx]
;    guard = guard[idx]
    
    if n_elements(vrng) eq 2 then begin
        idx = where(vs ge min(vrng) and vs le max(vrng))
        vs = vs[idx]
        ibias = ibias[idx]
    endif

    
    ofn = shomedir()+'/bias_iv_rbsp'+rbx+'_probe'+pair+'_'+time_string(tr[0],tformat='YYYY_MMDD')+'.pdf'
    sgopen, ofn, xsize = 6, ysize = 4, /inch
    plot, vs, ibias, $
        xrange = vrng, xstyle = 1, $
        psym = 1, symsize = 0.5, xtitle = '(V'+strmid(pair,0,1)+'+V'+strmid(pair,1,1)+')/2 (V)', ytitle = 'Ibias (nA)', $
        title = 'RBSP-'+strupcase(rbx)+', probe'+pair+', I-V, '+time_string(tr[0],tformat='YYYY_MMDD')
    sgclose
end


tr = time_double(['2014-04-29/00:00','2014-04-30/00:00'])
tr = time_double(['2014-04-18/06:00','2014-04-18/10:00'])

tr = time_double(['2015-11-14/00:00','2015-11-14/09:00'])
vrng = [-10,0]
for i = 0, 2 do begin
    test_vi, tr, 'a', 2*i+1, vrange = vrng
    test_vi, tr, 'b', 2*i+1, vrange = vrng
endfor

end