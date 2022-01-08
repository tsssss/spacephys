
pro bias_param_info, tvar, tr, info

    get_data, tvar, t0, d0
    if n_elements(tr) eq 0 then message, 'no time info ...'
    
    ; analyze data in given time range.
    idx = where(t0 ge tr[0] and t0 le tr[1])
    t1 = t0[idx]
    d1 = d0[idx]
    
;    idx = uniq(t1)
;    t1 = t1[idx]
;    d1 = d1[idx]
;    nrec = n_elements(d1)
    
;    idx = where(abs(d1-smode(d1)) gt stddev(d1))
;    t1 = t1[idx]
;    d1 = d1[idx]
    nrec = n_elements(d1)
    
    ; locate where value changes are.
    diff = abs(d1[1:nrec-1]-d1[0:nrec-2])
    avg = abs(mean(diff))
    std = stddev(diff)
    case 1 of
        std/avg gt 5: idx = where(diff ge std, cnt)
        std/avg gt 1: idx = where(diff gt avg, cnt)
        else: idx = where(diff ge min(diff), cnt)
    endcase
    if idx[0] ne 0 then idx = [idx[0]-1,idx]
    if idx[-1] ne nrec then idx = [idx,idx[-1]+1]
    nparam = n_elements(idx)-1
    if nparam le 0 then message, 'no paramemter found ...'
    ts = t1[idx]   ; times for change.
    vs = dblarr(nparam)
    for i = 0, nparam-1 do vs[i] = median(d1[idx[i]:idx[i+1]-1])
        
    info = {time:ts,value:vs,datarate:sdatarate(ts)}
;    print, info.value
    
;    store_data, tvar+'_info', ts, vs, info
end


pro test_bias, date, rbx, probe, reload = reload

    ; date. 'yyyy-mm-dd'
    ; rbx. 'a' or 'b'.
    ; probe. 1,2,3,4,5,6.

    dt = 86400d     ; sec.
    rbspx = 'rbsp'+rbx
    prefix = rbspx+'_'
    
    ; probes, can be [1,2], [3,4], [5,6].
    suffix = string(probe,format='(I0)')
    
    ; rough time range that includes the bias sweep for ALL probes.
    ; bias sweep for rbsp is always scheduled around the first apogee of the day.
    ; therefore, +/- 0.5 day, or 12 hr is enough to include all possible apogee.
    tr = time_double(date,tformat='YYYY-MM-DD')+dt*0.5*[-1,1]
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

            tvar = prefix+'E'+string(22*i+12,format='(I2)')
            store_data, tvar, t0, 0.5*(vsvy[*,2*i+1]-vsvy[*,2*i])*10
            options, tvar, 'ysubtitle', '[mV/m]'
            options, tvar, 'colors', i*2+2
            options, tvar, 'labels', 'E'+string(22*i+12,format='(I2)')

        endfor
    endif
    
    ; 'rbspx_IEFI_[IBIAS,GAURD,USHER][123456]' save the original parameters.
    ; 'rbspx_[ibias,guard,usher][123456]' save the stepped version.

    ; use usher and guard to locate rought time range.
    tvar = prefix+'guard'+suffix
    get_data, tvar, t0, guard
    tvar = prefix+'usher'+suffix
    get_data, tvar, t0, usher
    err = 0.5   ; V.
    idx = where(abs(guard-smode(guard)) gt err and abs(usher-smode(usher)) gt err)
    ttr = minmax(t0[idx])
    
    ; must have changing usher within the non-background times.
    diff = idx[1:*]-idx[0:-2]
    seps = where(diff gt 100, cnt) ; allow skip up to 100 records.
    if seps[0] ne 0 then seps = [0,seps]
    if seps[cnt-1] ne n_elements(idx)-1 then seps = [seps,n_elements(idx)-1]
    nsep = n_elements(seps)-1
    flags = bytarr(nsep)
    for i = 0, nsep-1 do begin
        idx1 = idx[seps[i]+1]
        idx2 = idx[seps[i+1]]
        diff = guard[idx1+1:idx2]-guard[idx1:idx2-1]
        if max(diff) le err then flags[i] = 1
    endfor
    tmp = where(flags eq 0, cnt)
    if cnt eq 0 then message, 'no sweep time found ...'
    idx1 = idx[seps[tmp[0]]+1]
    idx2 = idx[seps[tmp[0]+1]]
    ttr = t0[[idx1,idx2]]

    tvar = prefix+'IEFI_IBIAS'+suffix
    get_data, tvar, t0
    idx = where(t0 gt min(ttr) and t0 le max(ttr), cnt)
    t0 = t0[idx]
    tmp = t0[1:cnt-1]-t0[0:cnt-2]
    err = 60 ; sec.
    dr = smode(tmp, error = err)
    idx = where(abs(tmp-dr) le err, cnt)
    t0 = t0[idx]
    ttr = minmax(t0)
    ttr[1]+= dr
    nrec = cnt
    ttr+= err*[-1,1]

    ; read voltage-current relation.
    vinfo = {v:0d,ibias:0d,guard:0d,usher:0d}
    vinfo = replicate(vinfo, nrec)
    
    tvar = prefix+'V'+suffix
    get_data, tvar, t1, vs
    
    tvar = prefix+'ibias'+suffix
    get_data, tvar, t2, ibias
    tvar = prefix+'usher'+suffix
    get_data, tvar, t2, usher
    tvar = prefix+'guard'+suffix
    get_data, tvar, t2, guard
    
    
    for i = 0, nrec-1 do begin
        idx = where(t1 ge t0[i] and t1 lt t0[i]+dr)
        vinfo[i].v = mean(vs[idx])
        tmp = min(t2-t0[i], /absolute, idx)
        vinfo[i].ibias = ibias[idx]
        vinfo[i].usher = usher[idx]
        vinfo[i].guard = guard[idx]
    endfor
    tvar = prefix+'IV'+suffix
    store_data, tvar, ttr, vinfo
    
    
    tvar = prefix+'Vsc'
    vars = prefix+['Vsc'+['12','34','56'],'V'+suffix]
    stplot_merge, vars, newname = tvar
    vars = prefix+'Vsc'+['12','34']
    vars = vars[where(strpos(vars,suffix) eq -1)]
    stplot_minmax, vars, tmp, /get, trange = ttr
    options, tvar, 'yrange', tmp
    options, tvar, 'colors', [6,4,2,0]
    options, tvar, 'labels', ['Vsc'+['12','34','56'],'V'+suffix]
    options, tvar, 'ysubtitle', '[V]'


    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', -1


    ; prepare V and E.
    tvar = prefix+'V_comb'
    vars = ['1','2','3','4',suffix]
    vars = vars[uniq(vars,sort(vars))]
    vars = prefix+'V'+vars
    stplot_merge, vars, newname = tvar, colors = indgen(n_elements(vars))
    options, tvar, 'yrange', [-10,0]
    options, tvar, 'ystyle', 1
    options, tvar, 'V!C[V]'
    
    tvar = prefix+'E_comb'
    vars = prefix+['E'+['12','34','56']]
    stplot_merge, vars, newname = tvar
    options, tvar, 'yrange', [-10,10]
    options, tvar, 'ystyle', 1
    
    vars = prefix+['V_comb','E_comb',['ibias','guard','usher']+suffix]
    tplot, vars, trange = ttr+err*[-1,1]
end


date = '2013-09-11'
date = '2014-12-03'
date = '2014-12-17'

rbx = 'b'
probe = 1

test_bias, date, rbx, probe;, /reload
end
